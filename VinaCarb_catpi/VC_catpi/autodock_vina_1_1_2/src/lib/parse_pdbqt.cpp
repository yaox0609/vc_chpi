/*

   Copyright (c) 2006-2010, The Scripps Research Institute
   Copyright (c) 2015, The University of Georgia

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute
           
   Modifications for Vina-Carb 1.0 By: Anita K. Nivedha <nivedha@uga.edu>
                                       The Woods' Lab
                                       Complex Carbohydrate Research Center
                                       The University of Georgia

*/

#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include "parse_pdbqt.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "parse_error.h"
#include "glylib.h"


std::vector<int> branch_atom1;
std::vector<int> branch_atom2;
std::vector<size_t*> glyco_info;
std::vector< std::vector<size_t*> > ligand_glyco_info;
int hetatm=0;
std::ofstream VC_log;

struct stream_parse_error {
	unsigned line;
	std::string reason;
	stream_parse_error(unsigned line_, const std::string& reason_) : line(line_), reason(reason_) {}
	parse_error to_parse_error(const path& name) const {
		return parse_error(name, line, reason);
	}
};


std::vector<parsed_atom> ligand_info;

void add_context(context& c, std::string& str) {
	c.push_back(parsed_line(str, boost::optional<sz>()));
}

std::string omit_whitespace(const std::string& str, sz i, sz j) {
	if(i < 1) i = 1;
	if(j < i-1) j = i-1; // i >= 1
	if(j < str.size()) j = str.size();

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	// omit trailing whitespace
	while(i <= j && std::isspace(str[j-1]))
		--j;

	VINA_CHECK(i-1 < str.size());
	VINA_CHECK(j-i+1 < str.size());

	return str.substr(i-1, j-i+1);
}

struct atom_syntax_error {
	std::string nature;
	atom_syntax_error(const std::string& nature_) : nature(nature_) {}
};

template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
	VINA_CHECK(i >= 1);
	VINA_CHECK(i <= j+1);
	if(j > str.size()) throw atom_syntax_error("The line is too short");

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))
		++i;

	const std::string substr = str.substr(i-1, j-i+1);
	try {
		return boost::lexical_cast<T>(substr);
	}
	catch(...) {
		throw atom_syntax_error(std::string("\"") + substr + "\" is not a valid " + dest_nature);
	}
}

parsed_atom parse_pdbqt_atom_string(const std::string& str){ 
	unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "atom number");
	vec coords(checked_convert_substring<fl>(str, 31, 38, "coordinate"),
			   checked_convert_substring<fl>(str, 39, 46, "coordinate"),
			   checked_convert_substring<fl>(str, 47, 54, "coordinate"));
	fl charge = 0;
	if(!substring_is_blank(str, 69, 76))
		charge = checked_convert_substring<fl>(str, 69, 76, "charge");
	std::string name = omit_whitespace(str, 78, 79);
	std::string resname;
	std::string atomname;
	std::string resnum;
	std::string chainID;
	sz ad = string_to_ad_type(name);
	resname=str.substr(17,3); 
	resnum=str.substr(23,4); 
	atomname=str.substr(12,4); 
	chainID=str.substr(21,1);
	parsed_atom tmp(ad, charge, resname, resnum, atomname, coords, number, chainID);
	if(is_non_ad_metal_name(name))
		tmp.xs = XS_TYPE_Met_D;
	if(tmp.acceptable_type())
		return tmp;
	else 
		throw atom_syntax_error(std::string("\"") + name + "\" is not a valid AutoDock type. Note that AutoDock atom types are case-sensitive.");
}

struct atom_reference {
	sz index;
	bool inflex;
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct movable_atom : public atom_vc {
	vec relative_coords;
	std::string resname;
        std::string resnum;
        std::string atomname;
        std::string chainID;
	movable_atom(const atom_vc& a, const vec& relative_coords_, sz& ad_, std::string resname_, std::string resnum_, std::string atomname_, std::string chainID_) : atom_vc(a) {
		relative_coords = relative_coords_;
		this->ad = ad_;
		this->resname = resname_;
		this->resnum = resnum_;
		this->atomname = atomname_;
		this->chainID = chainID_;
	}
};

struct rigid {
	atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
	vector_mutable<ligand> ligands;
	vector_mutable<residue_vc> flex;

	mav atoms;
	atomv inflex;
	std::vector<parsed_atom> all_atoms_easy_access;  //Yao added 20230610

	distance_type_matrix atoms_atoms_bonds;
	matrix<distance_type> atoms_inflex_bonds;
	distance_type_matrix inflex_inflex_bonds;

	distance_type_matrix mobility_matrix() const {
		distance_type_matrix tmp(atoms_atoms_bonds);
		tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
		return tmp;
	}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
	struct node_t {
		sz context_index;
		parsed_atom a;
		std::vector<T> ps;
		node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

		// inflex atom insertion
		void insert_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
			nr.inflex.push_back(a);
		}
		void insert_immobiles_inflex(non_rigid_parsed& nr) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile_inflex(nr);
		}

		// insertion into non_rigid_parsed
		void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
			vec relative_coords; relative_coords = a.coords - frame_origin;
			c[context_index].second = nr.atoms.size();
			nr.atoms.push_back(movable_atom(a, relative_coords, a.ad, a.resname, a.resnum, a.atomname, a.chainID));
		}
		void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
			VINA_FOR_IN(i, ps)
				ps[i].insert_immobile(nr, c, frame_origin);
		}
	};

	typedef node_t<parsing_struct> node; 
	boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
	boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;  
	std::vector<parsed_atom> easy_access_atoms;

	void add(const parsed_atom& a, const context& c) { 
		VINA_CHECK(c.size() > 0);
		atoms.push_back(node(a, c.size()-1)); 
		easy_access_atoms.push_back(a);
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr);
		}
	}

	// insertion into non_rigid_parsed
	void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.atoms.size(), false);
			atoms[immobile_atom.get()].insert(nr, c, frame_origin);
		}
	}

	bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
		VINA_FOR_IN(i, atoms) {
			if(immobile_atom && immobile_atom.get() != i)
				return false;
			const node& nd = atoms[i];
			if(!nd.ps.empty())
				return false; // FIXME : iffy
		}
		return true;
	}
};

int BFMP(coord_3D **ring)
{
int RING=0, valid_RING_1C4_count=0,valid_RING_4C1_count=0;
coord_3D **four;
int i=0,j=0,k=0,no=4,cut_off=0,cut_off1=0,cut_check=0, echeck=0, fsize=0, no_planes=0, bpsize=0, sdsize=0, *bestplanes, *sortedplanes, *fiveatom_planes, sixcheck=0, fiveplanes=0, fivecheck2=0;
const char *conformers[16] = {"q","t","q","d","t","q","t","d","t","q","q","t","d","t","q"};
int list[60] = {0,1,2,3,0,1,2,4,0,1,2,5,0,1,3,4,0,1,3,5,0,1,4,5,0,2,3,4,0,2,3,5,0,2,4,5,0,3,4,5,1,2,3,4,1,2,3,5,1,2,4,5,1,3,4,5,2,3,4,5}; /*60 because when you have 4 atoms making up a dihedral, there are 4 ways in which those 4 atoms can be used to calculate the torsion sequentially*/
int six_list[25]= {0,1,2,3,1,2,3,4,2,3,4,5,3,4,5,0,4,5,0,1,5,0,1,2};
int secondlist[30] =   {4,5,3,5,3,4,2,5,2,4,2,3,1,5,1,4,1,3,1,2,0,5,0,4,0,3,0,2,0,1};

char *path_can, *path_vc, plane_output[4];
double dihedral, *fifteen_dihedrals, *d, *tendihedrals, *sortdihedrals, temp_dihedral,temp_dihedral1,temp_dihedral2,temp_dihedral3,  *six_dihedrals;
plane *fifteen_planes;
fileset Caset;

path_can=(char*)calloc(1000,sizeof(char));
bestplanes=(int*)calloc(0,sizeof(int));
sortedplanes=(int*)calloc(0,sizeof(int));
tendihedrals=(double*)calloc(0,sizeof(double));
sortdihedrals=(double*)calloc(0,sizeof(double));
fiveatom_planes=(int*)calloc(0,sizeof(int));
d=(double*)calloc(6,sizeof(double));
fifteen_dihedrals=(double*)calloc(16,sizeof(double));
fifteen_planes=(plane*)calloc(16,sizeof(plane));
four=(coord_3D**)calloc(60,sizeof(coord_3D*));
six_dihedrals=(double*)calloc(7,sizeof(double));

path_vc=getenv("VINA_CARB");
sprintf(path_can,"%s/autodock_vina_1_1_2/src/lib/canonicals.txt",path_vc);

if(cut_check==0){
        cut_off=10; //set default cutoff
        cut_off1=0-cut_off;
}

        for(i=0;i<15;i++){ //15 dihedrals, because 6 total atoms form a ring and 4 atoms are required to form a dihedral and 6c4 combinations equalts 360/24 = 15
                temp_dihedral=get_dihedral_ABCD_points(ring[list[j]][0],ring[list[j+1]][0],ring[list[j+2]][0],ring[list[j+3]][0]);
                temp_dihedral1=get_dihedral_ABCD_points(ring[list[j+1]][0],ring[list[j+2]][0],ring[list[j+3]][0],ring[list[j]][0]);
                temp_dihedral2=get_dihedral_ABCD_points(ring[list[j+2]][0],ring[list[j+3]][0],ring[list[j]][0],ring[list[j+1]][0]);
                temp_dihedral3= get_dihedral_ABCD_points(ring[list[j+3]][0],ring[list[j]][0],ring[list[j+1]][0],ring[list[j+2]][0]);
                fifteen_dihedrals[i]=(fabs(temp_dihedral)+fabs(temp_dihedral1)+fabs(temp_dihedral2)+fabs(temp_dihedral3))/4*(180/PI); //the average of the 4 combinations
                j=j+4;

        }
        i=0;
        j=0;



        for(i=0;i<15;i++){
                dihedral=get_dihedral_ABCD_points(ring[list[j]][0],ring[list[j+1]][0],ring[list[j+2]][0],ring[list[j+3]][0]);
                four[0]=&ring[list[j]][0];
                four[1]=&ring[list[j+1]][0];
                four[2]=&ring[list[j+2]][0];
                four[3]=&ring[list[j+3]][0];
                fifteen_planes[i]=get_plane_for_ring(no,four);
                j=j+4;
        }

        /*filtering all the planes with dihedral angles less than 10*/

        for(i=0;i<15;i++){
                if(fifteen_dihedrals[i]<=cut_off && fifteen_dihedrals[i]>=cut_off1){
                        bpsize++;
                        bestplanes=(int*)realloc(bestplanes,(bpsize*sizeof(int)));
                        bestplanes[bpsize-1]=i;
                        tendihedrals=(double*)realloc(tendihedrals,(bpsize*sizeof(double)));
                        tendihedrals[bpsize-1]=fifteen_dihedrals[i];
                        sortdihedrals=(double*)realloc(sortdihedrals,(bpsize*sizeof(double)));
                        sortdihedrals[bpsize-1]=fifteen_dihedrals[i];
                        no_planes++;
                }
        }

        double temp;
        for(i=0;i<bpsize;i++){
                for(j=i;j<bpsize;j++){
                        if(fabs(sortdihedrals[i]) > fabs(sortdihedrals[j])){
                                temp=sortdihedrals[i];
                                sortdihedrals[i]=sortdihedrals[j];
                                sortdihedrals[j]=temp;
                        }
                }
        }

        for(i=0;i<bpsize;i++){
                for(j=0;j<bpsize;j++){
                        if(sortdihedrals[i]==tendihedrals[j]){
                                sdsize++;
                                sortedplanes=(int*)realloc(sortedplanes,(sdsize*sizeof(int)));
                                sortedplanes[sdsize-1]=bestplanes[j];
                        }
                }
        }
        /*calculating the six dihedrals around the ring to check if it is a flat ring or an envelope*/

        j=0;
        for(i=0;i<6;i++){
                dihedral=get_dihedral_ABCD_points(ring[six_list[j]][0],ring[six_list[j+1]][0],ring[six_list[j+2]][0],ring[six_list[j+3]][0]);
                six_dihedrals[i]=dihedral*(180/PI);
                if(six_dihedrals[i]<=5 && six_dihedrals[i]>=-5){
                        sixcheck++;
                }
                j=j+4;
        }


int envelope_check=0;
int pre_e_check=0;
        if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[3])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[2])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[4])-fabs(six_dihedrals[5])) <= 9.0){
                pre_e_check++;
        }

if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[4])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[3])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[5])-fabs(six_dihedrals[0])) <= 9.0){
                pre_e_check++;
        }
}

if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[5])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[3])-fabs(six_dihedrals[4])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[1])) <= 9.0){
                pre_e_check++;
        }
}


if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[3])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[4])-fabs(six_dihedrals[5])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[2])) <= 9.0){
                pre_e_check++;
	}
}
if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[4])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[5])-fabs(six_dihedrals[0])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[3])) <= 9.0){
                pre_e_check++;
        }
}

if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[5])) <= 6.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[1])) <= 5.0){
                pre_e_check++;
        }
         if(fabs(fabs(six_dihedrals[3])-fabs(six_dihedrals[4])) <= 9.0){
                pre_e_check++;
        }
}

if(pre_e_check==3){
        envelope_check=1;
}
        int *envelope_planes;
        envelope_planes=(int*)calloc(0,sizeof(int));
        int ssize;
        int esize=0;
fiveplanes=0;
fivecheck2=0;
if(envelope_check==1){
         for(i=0;i<1;i++){
                        if((six_dihedrals[0] <=12.0 && six_dihedrals[0] >=-12.0) && (six_dihedrals[1]<=12.0 && six_dihedrals[1]>=-12.0)){
                                fivecheck2=1;
                                ssize=0;
                                for(j=ssize;j<ssize+8;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                        if((six_dihedrals[1] <=12.0 && six_dihedrals[1] >=-12.0) && (six_dihedrals[2]<=12.0 && six_dihedrals[2]>=-12.0)){
                                fivecheck2=1;
                                ssize=4;
                                for(j=ssize;j<ssize+8;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                        if((six_dihedrals[2] <=12.0 && six_dihedrals[2] >=-12.0) && (six_dihedrals[3]<=12.0 && six_dihedrals[3]>=-12.0)){
                                fivecheck2=1;
                                ssize=8;
                                for(j=ssize;j<ssize+8;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                        if((six_dihedrals[3] <=12.0 && six_dihedrals[3] >=-12.0) && (six_dihedrals[4]<=12.0 && six_dihedrals[4]>=-12.0)){
                                fivecheck2=1;
                                ssize=12;
                                for(j=ssize;j<ssize+8;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                        if((six_dihedrals[4] <=12.0 && six_dihedrals[4] >=-12.0) && (six_dihedrals[5]<=12.0 && six_dihedrals[5]>=-12.0)){
                                fivecheck2=1;
                                ssize=16;
                                for(j=ssize;j<ssize+8;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                        if((six_dihedrals[5] <=12.0 && six_dihedrals[5] >=-12.0) && (six_dihedrals[0]<=12.0 && six_dihedrals[0]>=-12.0)){
                                fivecheck2=1;
                                ssize=20;
                                for(j=ssize;j<ssize+4;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                ssize=0;
                                for(j=ssize;j<ssize+4;j++){
                                        esize++;
                                        envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                                        envelope_planes[esize-1]=six_list[j];

                                }
                                printf("yes\n");
                                fiveplanes++;
                        }
                }

}
i=0;
j=0;
if(envelope_check==1 && fivecheck2==1){
        while(i<fiveplanes*2*4){
                j=i;
                while(j<i+4){
                        echeck=0;
                        for(k=i+4;k<i+8;k++){
                                if(envelope_planes[j] == envelope_planes[k]){
                                        echeck=1;
                                }
                        }
                        if(echeck==0){
                                fsize++;
                                fiveatom_planes=(int*)realloc(fiveatom_planes,(fsize*sizeof(fiveatom_planes)));
                                fiveatom_planes[fsize-1]=envelope_planes[j];
                                for(k=i+4;k<i+8;k++){
                                        fsize++;
                                        fiveatom_planes=(int*)realloc(fiveatom_planes,(fsize*sizeof(fiveatom_planes)));
                                        fiveatom_planes[fsize-1]=envelope_planes[k];
                                }
                                j=i+4;
                        }

                        j++;
                }

                i=i+8;
        }//while(i<fiveplanes*2*4)

        int ringatoms[6]= {0,1,2,3,4,5};
        ssize=0;
        int *sixth_atom;
        sixth_atom=(int*)calloc(0,sizeof(sixth_atom));
        i=0;
        while(i<fiveplanes*5){
                j=i;
                k=0;
                while(k<6){
                        echeck=0;
                        for(j=i;j<i+5;j++){
                                if(fiveatom_planes[j]==ringatoms[k]){
                                        echeck=1;
                                }
                        }

                        if(echeck==0){
                                ssize++;
                                sixth_atom=(int*)realloc(sixth_atom,(ssize*sizeof*(sixth_atom)));
                                sixth_atom[ssize-1]=ringatoms[k];

                        }
                        k++;
                }

                i=i+5;
        }//while (i<fiveplanes*5)
        coord_3D **five;
        five=(coord_3D**)calloc(6,sizeof(coord_3D*));
        plane *fiveatomplane;
        fiveatomplane=(plane*)calloc(1,sizeof(plane));
        int n=5;

        for(i=0;i<fiveplanes;i++){
                five[0]=&ring[fiveatom_planes[i]][0];
                five[1]=&ring[fiveatom_planes[i+1]][0];
                five[2]=&ring[fiveatom_planes[i+2]][0];
                five[3]=&ring[fiveatom_planes[i+3]][0];
                five[4]=&ring[fiveatom_planes[i+4]][0];
                fiveatomplane[0]=get_plane_for_ring(n,five);
                d[i]=get_signed_distance_from_point_to_plane(fiveatomplane[0],ring[sixth_atom[i]][0]);
        }
}

if(fivecheck2==0 && no_planes==0){
}
Caset.N=strdup(path_can);
fileslurp Cslurp;
Cslurp=slurp_file(Caset);
char cname[6];
int dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9,dum10,dum11,dum12;
int *atoms;
atoms=(int*)calloc(13,sizeof(int));
int match=0,match1=0,ccheck=0,bcheck=0;
if(fivecheck2==0 && no_planes>0){
        if(sdsize>=3){
        for(i=0;i<Cslurp.n;i++){
                sscanf(Cslurp.L[i],"%s",cname);
                if(strstr(cname,"chair") != NULL){
                        sscanf(Cslurp.L[i],"%s%d%d%d%d%d%d%d%d%d%d%d%d",cname,&dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&dum8,&dum9,&dum10,&dum11,&dum12);
                        atoms[0]=dum1;atoms[1]=dum2;atoms[2]=dum3;atoms[3]=dum4;atoms[4]=dum5;atoms[5]=dum6;atoms[6]=dum7;atoms[7]=dum8;atoms[8]=dum9;atoms[9]=dum10;atoms[10]=dum11;atoms[11]=dum12;
                }
        }//for  
        }
        while(j<3){
                k=0;
                while(k<12){
                        match=0;
                        if(list[sortedplanes[j]*4] == atoms[k]){
                                match++;
                        }
                        if(list[sortedplanes[j]*4+1]== atoms[k+1]){
                                match++;
                        }
                        if(list[sortedplanes[j]*4+2] == atoms[k+2]){
                                match++;
                        }
                        if(list[sortedplanes[j]*4+3] == atoms[k+3]){
                               match++;
                        }
                        if(match==4){
                                k=12;
                        }
                        k=k+4;

                }//while k
                if(match==4){
                        match1++;
                }
                j++;
        }//while no_planes
       if(match1==3){
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring[3][0]);
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring[0][0]);
                if(d[0]<0.0 && d[1] > 0.0){
                        RING=2;
                }
                if(d[0]>0.0 && d[1] <0.0){
                        RING=1;
                }

                ccheck=1;
                }//match1==3    
        if(ccheck==0 && sdsize >=3){
                match=0;
                match1=0;
                int bnum=1;
                j=0;k=0;
                for(i=0;i<Cslurp.n;i++){
                        char conf_name[6];
                        sprintf(conf_name,"boat%d",bnum);
                        sscanf(Cslurp.L[i],"%s",cname);
                        if(strstr(cname,conf_name) != NULL){
                                sscanf(Cslurp.L[i],"%s%d%d%d%d%d%d%d%d%d%d%d%d",cname,&dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&dum8,&dum9,&dum10,&dum11,&dum12);
                                atoms[0]=dum1;atoms[1]=dum2;atoms[2]=dum3;atoms[3]=dum4;atoms[4]=dum5;atoms[5]=dum6;atoms[6]=dum7;atoms[7]=dum8;atoms[8]=dum9;atoms[9]=dum10;atoms[10]=dum11;atoms[11]=dum12;
                                bnum++;
                                j=0;
                                match1=0;
                                while(j<sdsize){
                                        k=0;
                                        match=0;
                                        while(k<12){
                                                match=0;
                                                if(list[sortedplanes[j]*4] == atoms[k]){
                                                        match++;
                                                }
                                                if(list[sortedplanes[j]*4+1] == atoms[k+1]){
                                                        match++;
                                                }
                                                if(list[sortedplanes[j]*4+2] == atoms[k+2]){
                                                        match++;
                                                }
                                                if(list[sortedplanes[j]*4+3] == atoms[k+3]){
                                                        match++;
                                                }
                                                if(match==4){
                                                        k=12;
                                                }
                                                k=k+4;
                                                }
                                                j++;
                                                if(match==4){
                                                        match1++;
                                                }

                                }//j<no_planes  
                        }
                        if(match1==3){
                                if(strstr(cname,"boat1") != NULL){
                                        d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring[3][0]);
                                        d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring[0][0]);
                                        if(d[0]>0.0 && d[1]>0.0){
                                                bcheck=1;
                                        }
                                        if(d[0]<0.0 && d[1]<0.0){
                                                bcheck=1;
                                        }

                                }
                                if(strstr(cname,"boat2") != NULL){
                                        d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[3],ring[2][0]);
                                        d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[3],ring[5][0]);
                                        if(d[0]>0.0 && d[1]>0.0){
                                                bcheck=1;
                                        }
                                        if(d[0]<0.0 && d[1]<0.0){
                                                bcheck=1;
                                        }

                                }
                                if(strstr(cname,"boat3") != NULL){
                                        d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[7],ring[1][0]);
                                        d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[7],ring[4][0]);
                                        if(d[0]>0.0 && d[1]>0.0){
                                                bcheck=1;
                                        }
                                        if(d[0]<0.0 && d[1]<0.0){
                                                bcheck=1;
                                        }
                                }
                                        i=Cslurp.n;
                        }

                }//for                          
        }//if ccheck=0;

        int scheck=0;
        int snum=1;
        int nsize=0;
        match1=0;
        match=0;
        if(ccheck==0 && bcheck==0){
        for(i=0;i<Cslurp.n;i++){
                char conf_name[6];
                sprintf(conf_name,"Skew%d",snum);
                sscanf(Cslurp.L[i],"%s",cname);
                if(strstr(cname,conf_name) != NULL){
                        sscanf(Cslurp.L[i],"%s%d%d%d%d%d%d%d%d",cname,&dum1,&dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&dum8);
                        atoms[0]=dum1;atoms[1]=dum2;atoms[2]=dum3;atoms[3]=dum4;atoms[4]=dum5;atoms[5]=dum6;atoms[6]=dum7;atoms[7]=dum8;
                        snum++;
                        j=0;
                        match1=0;
                        if(sdsize>=3){
                                nsize=3;
                        }
                        else{
                                nsize=sdsize;
                        }
                        while(j<nsize){
                                k=0;
                                match=0;
                                while(k<8){
                                        match=0;
                                        if(list[sortedplanes[j]*4] == atoms[k]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+1] == atoms[k+1]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+2] == atoms[k+2]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+3] == atoms[k+3]){
                                                match++;
                                        }
                                        if(match==4){
                                                k=8;
                                        }
                                        k=k+4;
                                 }
                                 j++;
                                 if(match==4){
                                        match1++;
                                 }

                          }//j<no_planes  
                  }
                  if(match1==2){
                        if(strstr(cname,"Skew1") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[11],ring[0][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[11],ring[4][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                        scheck=1;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                        scheck=1;
                                }

                        }
                        if(strstr(cname,"Skew2") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[6],ring[5][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[6],ring[1][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                        scheck=1;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                }
                        }
                        if(strstr(cname,"Skew3") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[13],ring[2][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[13],ring[0][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                        scheck=1;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                        scheck=1;
                                }
                        }
                        i=Cslurp.n;
                   }

         }//for                    
        }
        int hnum=1;
        nsize=0;
        match1=0;
        match=0;
        if(ccheck==0 && bcheck==0 && scheck==0){
        for(i=0;i<Cslurp.n;i++){
                char conf_name[6];
                sprintf(conf_name,"Half%d",hnum);
                sscanf(Cslurp.L[i],"%s",cname);
                if(strstr(cname,conf_name) != NULL){
                        sscanf(Cslurp.L[i],"%s%d%d%d%d",cname,&dum1,&dum2,&dum3,&dum4);
                        atoms[0]=dum1;atoms[1]=dum2;atoms[2]=dum3;atoms[3]=dum4;
                        hnum++;
                        j=0;
                        match1=0;
                        nsize=sdsize;
                        while(j<1){
                                k=0;
                                match=0;
                                while(k<4){
                                        match=0;
                                        if(list[sortedplanes[j]*4] == atoms[k]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+1] == atoms[k+1]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+2] == atoms[k+2]){
                                                match++;
                                        }
                                        if(list[sortedplanes[j]*4+3] == atoms[k+3]){
                                                match++;
                                        }
                                        if(match==4){
                                                k=4;
                                        }
                                        k=k+4;
                                }
                                j++;
                                if(match==4){
                                        match1++;
                                }

                        }//j<no_planes  
              }
                if(match1==1){
                        if(strstr(cname,"Half1") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[14],ring[0][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[14],ring[1][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                        i=Cslurp.n;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                        i=Cslurp.n;
                                }
                        }
                        if(strstr(cname,"Half2") != NULL){
                                 d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[9],ring[2][0]);
                                 d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[9],ring[1][0]);
                                 if(d[0]<0.0 && d[1]>0.0){
                                        i=Cslurp.n;
                                 }
                                 if(d[0]>0.0 && d[1]<0.0){
                                        i=Cslurp.n;
                                  }
                        }
                        if(strstr(cname,"Half3") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[5],ring[2][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[5],ring[3][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                        i=Cslurp.n;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                        i=Cslurp.n;
                                }
                       }
                       if(strstr(cname,"Half4") != NULL){
                                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[2],ring[4][0]);
                                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[2],ring[3][0]);
                                if(d[0]<0.0 && d[1]>0.0){
                                          i=Cslurp.n;
                                }
                                if(d[0]>0.0 && d[1]<0.0){
                                          i=Cslurp.n;
                                }
                       }
                       if(strstr(cname,"Half5") != NULL){
                                 d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[0],ring[4][0]);
                                 d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[0],ring[5][0]);
                                 if(d[0]<0.0 && d[1]>0.0){
                                            i=Cslurp.n;
                                 }
                                 if(d[0]>0.0 && d[1]<0.0){
                                            i=Cslurp.n;
                                 }
                        }
                        if(strstr(cname,"Half6") != NULL){
                                  d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[10],ring[0][0]);
                                  d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[10],ring[5][0]);
                                  if(d[0]<0.0 && d[1]>0.0){
                                             i=Cslurp.n;
                                  }
                                  if(d[0]>0.0 && d[1]<0.0){
                                             i=Cslurp.n;
                                  }
                        }
                  }//if match==1
        }//for                       
        }
                for(i=0;i<no_planes;i++){
                        d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[sortedplanes[i]],ring[secondlist[sortedplanes[i]*2]][0]);
                        d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[sortedplanes[i]],ring[secondlist[sortedplanes[i]*2+1]][0]);
                        if(d[0]<0.0 && d[1]>0.0){
			sprintf(plane_output,"%d%s%d",secondlist[sortedplanes[i]*2+1]+1,conformers[sortedplanes[i]],secondlist[sortedplanes[i]*2]+1);
                        }
                        if(d[0]>0.0 && d[1]<0.0){
                        sprintf(plane_output,"%d%s%d",secondlist[sortedplanes[i]*2]+1,conformers[sortedplanes[i]],secondlist[sortedplanes[i]*2+1]+1);
                        }
                        if(d[0]>0.0 && d[1] >0.0){
                        sprintf(plane_output,"%d%d%s",secondlist[sortedplanes[i]*2]+1,secondlist[sortedplanes[i]*2+1]+1,conformers[sortedplanes[i]]);
                        }
                        if(d[0]<0.0 && d[1]<0.0){
                        sprintf(plane_output,"%s%d%d",conformers[sortedplanes[i]],secondlist[sortedplanes[i]*2]+1,secondlist[sortedplanes[i]*2+1]+1);
                        }
if(strstr(plane_output,"2d5")!='\0'||strstr(plane_output,"4d1")!='\0'||strstr(plane_output,"6d3")!='\0')
	{
valid_RING_4C1_count++;
	}
if(strstr(plane_output,"5d2")!='\0'||strstr(plane_output,"1d4")!='\0'||strstr(plane_output,"3d6")!='\0')
	{
valid_RING_1C4_count++;
	}
                }
plane_output[0]='\0';
plane_output[1]='\0';
plane_output[2]='\0';
plane_output[3]='\0';

}//fivecheck2==0
if(RING==0 && valid_RING_1C4_count>=2)
{
RING=2;
}
else if (RING==0 && valid_RING_4C1_count>=2)
{
RING=1;
}
return RING;
}

unsigned parse_one_unsigned(const std::string& str, const std::string& start, unsigned count) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp;
	in_str >> tmp;
	if(!in_str || tmp < 0) 
		throw stream_parse_error(count, "Syntax error");
	return unsigned(tmp);
}

void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned count, unsigned& first, unsigned& second) {
	std::istringstream in_str(str.substr(start.size()));
	int tmp1, tmp2;
	in_str >> tmp1;
	in_str >> tmp2;
	if(!in_str || tmp1 < 0 || tmp2 < 0) 
		throw stream_parse_error(count, "Syntax error");
	first = unsigned(tmp1);
	second = unsigned(tmp2);
}

void parse_pdbqt_rigid(const path& name, rigid& r) {
	ifile in(name);
	unsigned count = 0;
	std::string str;

	while(std::getline(in, str)) {
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "TER")) {} // ignore 
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				//r.atoms.push_back(parse_pdbqt_atom_string(str)); //This is the original source code. 
				parsed_atom parsed = parse_pdbqt_atom_string(str);
				r.atoms.push_back(parsed);
				explicit_atom_update(parsed, r.atoms.back());
			}
			catch(atom_syntax_error& e) {
				throw parse_error(name, count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw parse_error(name, count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
}

double get_torsion_coords_vec_list(vec A, vec B, vec C, vec D)
{ //returns the torsion angle formed by 4 x,y,z co-ordinates.
double angle=0.0;
vec AB, BC, CD, ABCcross, BCDcross;
double ABC_BCD_dot, AB_BCD_dot, BCscalar_AB_BCD_dot;
AB=B-A;
BC=C-B;
CD=D-C;
ABCcross=cross_product(AB,BC);
BCDcross=cross_product(BC,CD);
ABC_BCD_dot=dot_product(ABCcross,BCDcross);
AB_BCD_dot=dot_product(AB,BCDcross);
BCscalar_AB_BCD_dot=magnitude(BC)*AB_BCD_dot;                        
angle=atan2(BCscalar_AB_BCD_dot,ABC_BCD_dot)*57.2957795; //angle in degrees
return angle;           
}                       

void parse_pdbqt_root_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				if(hetatm==1){
				parsed_atom LIGpa=parse_pdbqt_atom_string(str);//AKN
				ligand_info.push_back(LIGpa); 
				}
				p.add(parse_pdbqt_atom_string(str), c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "ENDROOT")) return;
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_root(std::istream& in, unsigned& count, parsing_struct& p, context& c) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "ROOT")) {
			parse_pdbqt_root_aux(in, count, p, c);
			break;
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux(std::istream& in, unsigned& count, const std::string& str, parsing_struct& p, context& c) {
	unsigned first, second;
	parse_two_unsigneds(str, "BRANCH", count, first, second); 
	sz i = 0;
	if(hetatm==1){
	branch_atom1.push_back(first);
	branch_atom2.push_back(second);
	}
	for(; i < p.atoms.size(); ++i)
		{
		if(p.atoms[i].a.number == first) {
			p.atoms[i].ps.push_back(parsing_struct());
			parse_pdbqt_branch(in, count, p.atoms[i].ps.back(), c, first, second);
			//Add easy access atoms from all child parsing structs to the parent.
			std::vector<parsed_atom>& child_easy_acess_atoms = p.atoms[i].ps.back().easy_access_atoms;
			p.easy_access_atoms.insert(p.easy_access_atoms.end(), child_easy_acess_atoms.begin(), child_easy_acess_atoms.end());
			break;
		}
		}
	if(i == p.atoms.size())
		throw stream_parse_error(count, "No atom number " + boost::lexical_cast<std::string>(first) + " in this branch");
}

              

void parse_pdbqt_aux(std::istream& in, unsigned& count, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue_vc) {
	parse_pdbqt_root(in, count, p, c);

	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(!residue_vc && starts_with(str, "TORSDOF")) {
			if(torsdof) throw stream_parse_error(count, "TORSDOF can occur only once");
			torsdof = parse_one_unsigned(str, "TORSDOF", count);
		}
		else if(residue_vc && starts_with(str, "END_RES")) return; 
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
	int branch=0, i=0;
	VINA_FOR(i,ligand_info.size())	
	{//This loop helps to remove blank spaces from atomnames in case they have them so the actual names can be used for comparison instead of with trailing blank spaces.
		if(ligand_info[i].atomname.find(" ")<3)
		{
		ligand_info[i].atomname.replace(ligand_info[i].atomname.find(" "),1,"\0");
		}
	}
	VINA_FOR(h,branch_atom1.size())
	{
	int ring1_conf=0, ring2_conf=0;
	size_t *S1_C1, *S1_O5, *S1_O1, *S1_C5, *S1_C2, *S1_C3, *S1_C4, *S2_Cxp1, *S2_C2, *S2_O5, *S2_C5, *S2_C1, *S2_O1, *S2_Ox, *S2_Cx, *S2_Cxm1, *S2_O4, *S2_C4, *S2_C3, *sizet_ring1_conf;
	size_t *S1_AB, *S2_AE, *S2_Link, *S2_6AE, *sizet_ring2_conf;
	coord_3D *C1, *C2, *C3, *C4, *C5, *O5, **ring_atoms;
        double S1_AB_angle=0.0, S2_AE_angle=0.0, S2_omega_angle=0.0;
	ring_atoms=(coord_3D**)calloc(6,sizeof(coord_3D*));
	C1=(coord_3D*)calloc(1,sizeof(coord_3D));
	C2=(coord_3D*)calloc(1,sizeof(coord_3D));
	C3=(coord_3D*)calloc(1,sizeof(coord_3D));
	C4=(coord_3D*)calloc(1,sizeof(coord_3D));
	C5=(coord_3D*)calloc(1,sizeof(coord_3D));
	O5=(coord_3D*)calloc(1,sizeof(coord_3D));
        S1_C1=(size_t*)calloc(1,sizeof(size_t));
        S1_O5=(size_t*)calloc(1,sizeof(size_t));
        S1_O1=(size_t*)calloc(1,sizeof(size_t));
        S1_C5=(size_t*)calloc(1,sizeof(size_t));
        S1_C2=(size_t*)calloc(1,sizeof(size_t));
        S1_C3=(size_t*)calloc(1,sizeof(size_t));
        S1_C4=(size_t*)calloc(1,sizeof(size_t));
        sizet_ring1_conf=(size_t*)calloc(1,sizeof(size_t));
        S2_Cxp1=(size_t*)calloc(1,sizeof(size_t));
        S2_C2=(size_t*)calloc(1,sizeof(size_t));
        S2_O5=(size_t*)calloc(1,sizeof(size_t));
        S2_O4=(size_t*)calloc(1,sizeof(size_t));
        S2_C4=(size_t*)calloc(1,sizeof(size_t));
        S2_C3=(size_t*)calloc(1,sizeof(size_t));
        S2_C5=(size_t*)calloc(1,sizeof(size_t));
        S2_C1=(size_t*)calloc(1,sizeof(size_t));
        S2_O1=(size_t*)calloc(1,sizeof(size_t));
        S2_Ox=(size_t*)calloc(1,sizeof(size_t));
        S2_Cx=(size_t*)calloc(1,sizeof(size_t));
        S2_Cxm1=(size_t*)calloc(1,sizeof(size_t));
        S1_AB=(size_t*)calloc(1,sizeof(size_t));
        S2_AE=(size_t*)calloc(1,sizeof(size_t));
        S2_Link=(size_t*)calloc(1,sizeof(size_t));
        S2_6AE=(size_t*)calloc(1,sizeof(size_t));
        sizet_ring2_conf=(size_t*)calloc(1,sizeof(size_t));
	S1_AB[0]=-1;
	S2_6AE[0]=-1;
        if((ligand_info[branch_atom1[h]-1].resnum!=ligand_info[branch_atom2[h]-1].resnum) && (ligand_info[branch_atom1[h]-1].resname.compare("OME")!=0 && ligand_info[branch_atom1[h]-1].resname.compare("ROH")!=0 ) && (ligand_info[branch_atom2[h]-1].resname.compare("OME")!=0 && ligand_info[branch_atom2[h]-1].resname.compare("ROH")!=0 )  )
	{//found glycosidic
	VINA_FOR(i,2)
		{//using 2 values for i, one for branch_atom1 and another for branch_atom2
			if(i==0)
			{
			branch=branch_atom1[h];
			}
			else if(i==1)
			{
			branch=branch_atom2[h];
			}
			if(ligand_info[branch-1].atomname.compare("C1")==0)
			{
			S1_C1[0]=branch-1;
				VINA_FOR(j,ligand_info.size())
				{
					if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
					{
						if(ligand_info[j].atomname.compare("C2")==0)
						{
						S1_C2[0]=j;
						}
						if(ligand_info[j].atomname.compare("O5")==0)
						{
						S1_O5[0]=j;
						}
						if(ligand_info[j].atomname.compare("C5")==0)
						{
						S1_C5[0]=j;
						}
						if(ligand_info[j].atomname.compare("C4")==0)
						{
						S1_C4[0]=j;
						}
						if(ligand_info[j].atomname.compare("C3")==0)
						{
						S1_C3[0]=j;
						}
					}
				}//end of j lig info size checking for C1
			}//end of if where C1 was found
			else if(ligand_info[branch-1].atomname[0]=='O')
			{
			S2_Ox[0]=branch-1;
				VINA_FOR(j,ligand_info.size())
                                {//in this loop I will find Cx, C2, C5 and C1
					if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
					{
						if(ligand_info[j].atomname[0]=='C')
						{
							if(ligand_info[j].atomname[1]==ligand_info[branch-1].atomname[1]) //found Cx
							{
							S2_Cx[0]=j;
								if(ligand_info[j].atomname[1]=='2')
								{//Finding the linkage type
								S2_Link[0]=2;
								}
								if(ligand_info[j].atomname[1]=='3')
								{
								S2_Link[0]=3;
								}
								if(ligand_info[j].atomname[1]=='4')
								{
								S2_Link[0]=4;
								}
								if(ligand_info[j].atomname[1]=='6')
								{
								S2_Link[0]=6;
								}
							}
							if(ligand_info[j].atomname.compare("C2")==0) //found Cx
							{
							S2_C2[0]=j;
							}
							if(ligand_info[j].atomname.compare("C5")==0) //found Cx
							{
							S2_C5[0]=j;
							}
							if(ligand_info[j].atomname.compare("C1")==0) //found Cx
							{
							S2_C1[0]=j;
							}
						}
					}
				} //end of j lig info size checking for Ox
			float distance=0.0, distance_1up=100.0, distance_1down=100.0;
				VINA_FOR(j,ligand_info.size())
                                {
				if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
                                        {
						if(ligand_info[j].atomname[0]=='C')
						{
						distance=sqrt(pow(ligand_info[j].coords.data[0]-ligand_info[S2_Cx[0]].coords.data[0],2)+pow(ligand_info[j].coords.data[1]-ligand_info[S2_Cx[0]].coords.data[1],2)+pow(ligand_info[j].coords.data[2]-ligand_info[S2_Cx[0]].coords.data[2],2));
							if(distance<distance_1up && ligand_info[j].atomname.compare(ligand_info[S2_Cx[0]].atomname)>0)
							{
							distance_1up=distance;
							S2_Cxp1[0]=j;
							}
							if(distance<distance_1down && ligand_info[j].atomname.compare(ligand_info[S2_Cx[0]].atomname)<0)
							{
							distance_1down=distance;
                                                        S2_Cxm1[0]=j;
							}
						}
					}
				}
				float distance_O1=100.0;
				VINA_FOR(j,ligand_info.size())
				{
					if(ligand_info[j].resnum==ligand_info[branch-1].resnum)
                                        {
						if(ligand_info[j].atomname.compare("O5")==0)
						{
						S2_O5[0]=j;
						}
						if(ligand_info[j].atomname.compare("O1")==0)
                                                {//To find S2_O1 -- if present inside same residue as O1
                                                S2_O1[0]=j;
                                                }
						if(ligand_info[j].atomname.compare("O4")==0)
                                                {//To find S2_O1 -- if present inside same residue as O1
                                                S2_O4[0]=j;
                                                }
                                                if(ligand_info[j].atomname.compare("C4")==0)
                                                {//To find S2_O1 -- if present inside same residue as O1
                                                S2_C4[0]=j;
                                                }
                                                if(ligand_info[j].atomname.compare("C3")==0)
                                                {//To find S2_O1 -- if present inside same residue as O1
                                                S2_C3[0]=j;
                                                }
						else
						{//To find S2_O1 -- when having to use Ox of neighbouring residue in it's place
							VINA_FOR(k,ligand_info.size())
							{//Going through all atoms im structure so that nearest linking Oxygen atom can be found
								if(ligand_info[k].atomname[0]=='O' && ligand_info[k].atomname[1]!='5')
								{
								distance=sqrt(pow(ligand_info[k].coords.data[0]-ligand_info[S2_C1[0]].coords.data[0],2)+pow(ligand_info[k].coords.data[1]-ligand_info[S2_C1[0]].coords.data[1],2)+pow(ligand_info[k].coords.data[2]-ligand_info[S2_C1[0]].coords.data[2],2));
									if(distance<distance_O1)
									{
									distance_O1=distance;
									S2_O1[0]=k;
									}
								}
							}
						}
					}
				}
			} //end of if where Ox was found
		}//end of if i is either 0 or 1!
				C1[0].i=ligand_info[S1_C1[0]].coords[0];
				C1[0].j=ligand_info[S1_C1[0]].coords[1];
				C1[0].k=ligand_info[S1_C1[0]].coords[2];

				C2[0].i=ligand_info[S1_C2[0]].coords[0];
				C2[0].j=ligand_info[S1_C2[0]].coords[1];
				C2[0].k=ligand_info[S1_C2[0]].coords[2];

				C3[0].i=ligand_info[S1_C3[0]].coords[0];
				C3[0].j=ligand_info[S1_C3[0]].coords[1];
				C3[0].k=ligand_info[S1_C3[0]].coords[2];

				C4[0].i=ligand_info[S1_C4[0]].coords[0];
				C4[0].j=ligand_info[S1_C4[0]].coords[1];
				C4[0].k=ligand_info[S1_C4[0]].coords[2];

				C5[0].i=ligand_info[S1_C5[0]].coords[0];
				C5[0].j=ligand_info[S1_C5[0]].coords[1];
				C5[0].k=ligand_info[S1_C5[0]].coords[2];

				O5[0].i=ligand_info[S1_O5[0]].coords[0];
				O5[0].j=ligand_info[S1_O5[0]].coords[1];
				O5[0].k=ligand_info[S1_O5[0]].coords[2];
				
                                S1_AB_angle=get_torsion_coords_vec_list(ligand_info[S1_C5[0]].coords,ligand_info[S1_O5[0]].coords,ligand_info[S1_C1[0]].coords,ligand_info[S2_Ox[0]].coords);
                                S2_omega_angle=get_torsion_coords_vec_list(ligand_info[S2_O4[0]].coords,ligand_info[S2_C4[0]].coords,ligand_info[S2_C3[0]].coords,ligand_info[S2_C5[0]].coords);
                                S2_AE_angle=get_angle_ABC(ligand_info[S2_Ox[0]].coords,ligand_info[S2_Cx[0]].coords,ligand_info[S2_O5[0]].coords);


				ring_atoms[0]=C1;
				ring_atoms[1]=C2;
				ring_atoms[2]=C3;
				ring_atoms[3]=C4;
				ring_atoms[4]=C5;
				ring_atoms[5]=O5;
				ring1_conf=BFMP(ring_atoms);
				sizet_ring1_conf[0]=ring1_conf;
				C1[0].i=ligand_info[S2_C1[0]].coords[0];
				C1[0].j=ligand_info[S2_C1[0]].coords[1];
				C1[0].k=ligand_info[S2_C1[0]].coords[2];

				C2[0].i=ligand_info[S2_C2[0]].coords[0];
				C2[0].j=ligand_info[S2_C2[0]].coords[1];
				C2[0].k=ligand_info[S2_C2[0]].coords[2];

				C3[0].i=ligand_info[S2_C3[0]].coords[0];
				C3[0].j=ligand_info[S2_C3[0]].coords[1];
				C3[0].k=ligand_info[S2_C3[0]].coords[2];

				C4[0].i=ligand_info[S2_C4[0]].coords[0];
				C4[0].j=ligand_info[S2_C4[0]].coords[1];
				C4[0].k=ligand_info[S2_C4[0]].coords[2];

				C5[0].i=ligand_info[S2_C5[0]].coords[0];
				C5[0].j=ligand_info[S2_C5[0]].coords[1];
				C5[0].k=ligand_info[S2_C5[0]].coords[2];

				O5[0].i=ligand_info[S2_O5[0]].coords[0];
				O5[0].j=ligand_info[S2_O5[0]].coords[1];
				O5[0].k=ligand_info[S2_O5[0]].coords[2];

                                ring_atoms[0]=C1;
                                ring_atoms[1]=C2;
                                ring_atoms[2]=C3;
                                ring_atoms[3]=C4;
                                ring_atoms[4]=C5;
                                ring_atoms[5]=O5;
				ring2_conf=BFMP(ring_atoms);
				sizet_ring2_conf[0]=ring2_conf;
				if(sizet_ring1_conf[0]==0)
				{
				VC_log<<"CHI energy penalties NOT applied to phi torsion in linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring2_conf[0]==0)
				{
				VC_log<<"CHI energy penalties NOT applied to psi torsion in linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring1_conf[0]==0 && sizet_ring2_conf[0]==0)
				{
				VC_log<<"CHI energy penalities NOT applied to glycosidic linkage b/w "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<".\n";
				}
				if(sizet_ring1_conf[0]!=0 && sizet_ring2_conf[0]!=0)
				{
				VC_log<<"CHI energy penalties applied to "<<ligand_info[branch_atom1[h]-1].resname<<" "<<ligand_info[branch_atom1[h]-1].resnum<<" and "<<ligand_info[branch_atom2[h]-1].resname<<" "<<ligand_info[branch_atom2[h]-1].resnum<<" linkage.\n";
				}
				if( ((S1_AB_angle<-40) && (S1_AB_angle>-80)) || ((S1_AB_angle>40) && (S1_AB_angle<80)) )
				{
					//Phi_Alpha_L
					S1_AB[0]=0;
				}
				else if( ((S1_AB_angle<-160) && (S1_AB_angle>-200)) || ((S1_AB_angle>160) && (S1_AB_angle<200)) )
				{
                                        //Phi_Beta_D
					S1_AB[0]=1;
				}
				if( (S2_AE_angle<120 && S2_AE_angle>80) )
				{
				//Axial attachment
				S2_AE[0]=0;
				}
				else if(S2_AE_angle>130 && S2_AE_angle<170)
				{
				//Equatorial attachment
				S2_AE[0]=1;
				}
                                if(S2_omega_angle>0)
                                        {
                                        S2_6AE[0]=1;
                                        }
                                else if (S2_omega_angle<0)
                                        {
                                        S2_6AE[0]=0;
                                        }
		//Co-ordinates START
		glyco_info.push_back(S1_O5); //index 0
		glyco_info.push_back(S1_C1); //index 1
		glyco_info.push_back(S2_Ox); //index 2
		glyco_info.push_back(S2_Cx); //index 3
		glyco_info.push_back(S2_Cxm1); //index 4
		//Co-ordinates END
		glyco_info.push_back(S1_AB); //index 5
		glyco_info.push_back(S2_AE); //index 6
		glyco_info.push_back(S2_Link); //index 7 
                glyco_info.push_back(S2_6AE); //index 8  //0 -> positive; 1-> negative
                //Co-ordinates START
                glyco_info.push_back(S2_O5); //index 9 
		glyco_info.push_back(sizet_ring1_conf); //index 10
		glyco_info.push_back(sizet_ring2_conf); //index 11
                //Co-ordinates ED
		ligand_glyco_info.push_back(glyco_info);

		glyco_info.clear();
		}//end of finding glycosidic
	}//end of for for branch atom size
}


std::vector<parsed_atom> liginfo_return()
{
return ligand_info;
}

std::vector< std::vector<size_t*> > glycan_info_func()
{
return ligand_glyco_info;
}

void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
	if(atm)
		VINA_RANGE(i, r.begin, r.end) {
			atom_reference& ar = atm.get();
			if(ar.inflex){ 
				nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
                        }
			else{
				nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
                       }
		}
}

void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
	if(axis_begin && axis_end) {
		atom_reference& r1 = axis_begin.get();
		atom_reference& r2 = axis_end  .get();
		if(r2.inflex) {
			VINA_CHECK(r1.inflex); // no atom-inflex rotors
			nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
		}
		else
			if(r1.inflex){
				nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
                        }
			else{
				nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
                            }
	}
}

typedef std::pair<sz, sz> axis_numbers;
typedef boost::optional<axis_numbers> axis_numbers_option;

void nr_update_matrixes(non_rigid_parsed& nr) {
	// atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)
	nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
	nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
	nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

template<typename B> // B == branch or main_branch or flexible_body 
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
	b.node.begin = nr.atoms.size();
	VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
		parsing_struct::node& p_node = p.atoms[i];
		if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
		else p_node.insert(nr, c, b.node.get_origin());
		p_node.insert_immobiles(nr, c, b.node.get_origin());
	}
	b.node.end = nr.atoms.size();
	nr_update_matrixes(nr);
	add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
	add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
	set_rotor(nr, p.axis_begin, p.axis_end);

	VINA_RANGE(i, b.node.begin, b.node.end){
		VINA_RANGE(j, i+1, b.node.end){
			nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME
                }
        }


	VINA_FOR_IN(i, p.atoms) { 	// postprocess children
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
				b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
				postprocess_branch(nr, ps, c, b.children.back());
			}
		}
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
	VINA_CHECK(!p.atoms.empty());
	nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
	postprocess_branch(nr, p, c, nr.ligands.back());
	nr_update_matrixes(nr); // FIXME ?
}

void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		p_node.insert_inflex(nr);
		p_node.insert_immobiles_inflex(nr);
	}
	VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
		parsing_struct::node& p_node = p.atoms[i];
		VINA_FOR_IN(j, p_node.ps) {
			parsing_struct& ps = p_node.ps[j];
			if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
				nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_//branch will assign begin and end
				postprocess_branch(nr, ps, c, nr.flex.back());
			}
		}
	}
	nr_update_matrixes(nr); // FIXME ?
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
	VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	VC_log.open("VC_log.txt");
	unsigned count = 0;
	parsing_struct p;
	boost::optional<unsigned> torsdof;
	try {
		parse_pdbqt_aux(in, count, p, c, torsdof, false);
		if(p.atoms.empty()) 
			throw parse_error(name, count, "No atoms in the ligand");
		if(!torsdof)
			throw parse_error(name, count, "Missing TORSDOF");
		postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint
	}
	catch(stream_parse_error& e) {
		throw e.to_parse_error(name);
	}

	nr.all_atoms_easy_access =  p.easy_access_atoms; //Yao added 20230610. Code is too complex, have no choice but doing this. 
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
	VC_log.close();
}

void parse_pdbqt_residue(std::istream& in, unsigned& count, parsing_struct& p, context& c) { 
	boost::optional<unsigned> dummy;
	parse_pdbqt_aux(in, count, p, c, dummy, true);
}

void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
	ifile in(name);
	unsigned count = 0;
	std::string str;

	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} // ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BEGIN_RES")) {
			try {
				parsing_struct p;
				parse_pdbqt_residue(in, count, p, c);
				postprocess_residue(nr, p, c);
			}
			catch(stream_parse_error& e) {
				throw e.to_parse_error(name);
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw parse_error(name, count, "Unknown or inappropriate tag");
	}
	VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_branch(std::istream& in, unsigned& count, parsing_struct& p, context& c, unsigned from, unsigned to) {
	std::string str;
	while(std::getline(in, str)) {
		add_context(c, str);
		++count;
		if(str.empty()) {} //ignore ""
		else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
		else if(starts_with(str, "REMARK")) {} // ignore
		else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, count, str, p, c);
		else if(starts_with(str, "ENDBRANCH")) {
			unsigned first, second;
			parse_two_unsigneds(str, "ENDBRANCH", count, first, second);
			if(first != from || second != to) 
				throw stream_parse_error(count, "Inconsistent branch numbers");
			if(!p.immobile_atom) 
				throw stream_parse_error(count, "Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch");
			return;
		}
		else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				parsed_atom a = parse_pdbqt_atom_string(str);
				
				if(hetatm==1){
				ligand_info.push_back(a);
				}
				if(a.number == to){
					p.immobile_atom = p.atoms.size();
                                }
				p.add(a, c);
			}
			catch(atom_syntax_error& e) {
				throw stream_parse_error(count, "ATOM syntax incorrect: " + e.nature);
			}
			catch(...) { 
				throw stream_parse_error(count, "ATOM syntax incorrect");
			}
		}
		else if(starts_with(str, "MODEL"))
			throw stream_parse_error(count, "Unexpected multi-MODEL input. Use \"vina_split\" first?");
		else throw stream_parse_error(count, "Unknown or inappropriate tag");
	}
}


//////////// new stuff //////////////////


struct pdbqt_initializer {
	model m;
	void initialize_from_rigid(const rigid& r) { // static really
		VINA_CHECK(m.grid_atoms.empty());
		m.grid_atoms = r.atoms;
		VINA_FOR_IN(i, m.grid_atoms){
			explicit_atom_update(r.atoms[i], m.grid_atoms[i]);
		}
	}
	void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
		VINA_CHECK(m.ligands.empty());
		VINA_CHECK(m.flex   .empty());

		m.ligands = nrp.ligands;
		m.flex    = nrp.flex;

		VINA_CHECK(m.atoms.empty());

		sz n = nrp.atoms.size() + nrp.inflex.size();
		m.atoms.reserve(n);
		m.coords.reserve(n);

		VINA_FOR_IN(i, nrp.atoms) {
			const movable_atom& a = nrp.atoms[i];
			atom_vc b = static_cast<atom_vc>(a);
			b.coords = a.relative_coords;
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
			explicit_atom_update(a, m.atoms.back());
		}
		VINA_FOR_IN(i, nrp.inflex) {
			const atom_vc& a = nrp.inflex[i];
			atom_vc b = a;
			b.coords = zero_vec_vc; // to avoid any confusion; presumably these will never be looked at
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
			explicit_atom_update(a, m.atoms.back());
		}

		VINA_CHECK(m.coords.size() == n);

		m.internal_coords.resize(m.coords.size(), zero_vec_vc); // FIXME

		m.minus_forces = m.coords;
		m.m_num_movable_atoms = nrp.atoms.size();

		if(is_ligand) {
			VINA_CHECK(m.ligands.size() == 1);
			m.ligands.front().cont = c;
		}
		else
			m.flex_context = c;

	}
	void initialize(const distance_type_matrix& mobility) {
		m.initialize(mobility);
	}
};

model parse_ligand_pdbqt  (const path& name) { // can throw parse_error
	non_rigid_parsed nrp;
	context c;
	hetatm=1;
	parse_pdbqt_ligand(name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_nrp(nrp, c, true);
	tmp.initialize(nrp.mobility_matrix());

	//Here, run something like tmp.m detect_rec_ar_rings. Yao added 20230602
	//tmp.m.build_residue_info();
	//tmp.m.build_lig_ar_ring_info();
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name, const path& flex_name) { // can throw parse_error
	hetatm=0;
	rigid r;
	non_rigid_parsed nrp;
	context c;
	parse_pdbqt_rigid(rigid_name, r);
	parse_pdbqt_flex(flex_name, nrp, c);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	tmp.initialize_from_nrp(nrp, c, false);
	tmp.initialize(nrp.mobility_matrix());

	//Here, run something like tmp.m.build_rec_ar_ring_info(). Yao added 20230602
	//tmp.m.build_residue_info();
	//tmp.m.build_rec_ar_ring_info();
	return tmp.m;
}

model parse_receptor_pdbqt(const path& rigid_name) { // can throw parse_error
	hetatm=0;
	rigid r;
	parse_pdbqt_rigid(rigid_name, r);

	pdbqt_initializer tmp;
	tmp.initialize_from_rigid(r);
	distance_type_matrix mobility_matrix;
	tmp.initialize(mobility_matrix);
	return tmp.m;
}
