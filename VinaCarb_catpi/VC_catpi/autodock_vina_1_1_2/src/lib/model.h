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

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <boost/optional.hpp> // for context
#include <map>
#include <algorithm> //std::find
#include <cmath>

#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "precalculate.h"
#include "igrid.h"
#include "grid_dim.h"

#include "atom.h"
struct interacting_pair {
	sz type_pair_index;
	sz a; 
	sz b;
	interacting_pair(sz type_pair_index_, sz a_, sz b_) : type_pair_index(type_pair_index_), a(a_), b(b_) {} 
};

typedef std::vector<interacting_pair> interacting_pairs; 

typedef std::pair<std::string, boost::optional<sz> > parsed_line; //Pair: This class couples together a pair of values, which may be of different types (T1 and T2). The individual values can be accessed through the public members first and second.
typedef std::vector<parsed_line> context; 

typedef std::vector<atom_vc*> aptrv; //Yao 20230604, vector of pointers to object atom_vc. Use pointers to avoid duplication of objects. 

struct ligand : public flexible_body, atom_range {
	unsigned degrees_of_freedom; // can be different from the apparent number of rotatable bonds, because of the disabled torsions
	interacting_pairs pairs;
	context cont;
	ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
	void set_range();
};

struct residue_vc : public main_branch {
	residue_vc(const main_branch& m) : main_branch(m) {}
};

//Yao added 20230610
struct model_residue {
	std::string name;
	std::string number;
	std::string chainID;
	aptrv atoms;
	std::map<atom_vc*, sz> aptr_index_map;
	model_residue(std::string& resname, std::string& resnum, std::string& chainID_){
		name = resname;
		number = resnum;
		chainID = chainID_;
	}
};
typedef std::vector<model_residue> resv;

//Yao added 20230706
typedef std::vector<vec*> vptrv;
struct ring_attribute {
	ring_attribute(szv& indices_, aptrv& atom_ptrs_){
		//VINA_CHECK(!this->atom_indices.empty());
		//VINA_CHECK(this->atom_indices.size() >= 3); //A ring must have >3 atoms.
		this->atom_indices = indices_;
		this->atom_ptrs = atom_ptrs_;
		this->atom_coord_ptrs.clear();
	}

	void ComputeCentroid(){
		fl avg_x = 0, avg_y = 0, avg_z = 0;
                VINA_FOR_IN(i, this->atom_coord_ptrs){
                        vec* coord = this->atom_coord_ptrs[i];

                        fl& x = coord->data[0]; fl& y = coord->data[1]; fl& z = coord->data[2];
                        avg_x += x; avg_y += y; avg_z += z;
                }

                fl rsize = (fl) this->atom_indices.size();
                avg_x /= rsize; avg_y /= rsize; avg_z /= rsize;
               	this->centroid.data[0] = avg_x; this->centroid.data[1] = avg_y; this->centroid.data[2] = avg_z;
		return;
	}

	void ComputeNormal(){
		vec *a1c = this->atom_coord_ptrs[0], *a2c = this->atom_coord_ptrs[1], *a3c = this->atom_coord_ptrs[2];
                fl V1x = a2c->data[0] - a1c->data[0], V1y = a2c->data[1] - a1c->data[1], V1z = a2c->data[2] - a1c->data[2];
                fl V2x = a3c->data[0] - a2c->data[0], V2y = a3c->data[1] - a2c->data[1], V2z = a3c->data[2] - a2c->data[2];
                this->normal.data[0] = V1y*V2z - V1z*V2y; this->normal.data[1] = V1z*V2x - V1x*V2z; this->normal.data[2] = V1x*V2y - V1y*V2x;
		normalize_vec_in_place(this->normal);
		return;
	}

	void ComputeEffectiveRadius(){
		fl total_distance = 0.00;
                VINA_FOR_IN(i, this->atom_coord_ptrs){
                        vec* coord = this->atom_coord_ptrs[i];
                        fl dx = this->centroid.data[0] - coord->data[0];
                        fl dy = this->centroid.data[1] - coord->data[1];
                        fl dz = this->centroid.data[2] - coord->data[2];
                        fl r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        total_distance  += r;
                }
                this->effective_radius =  total_distance / (fl) this->atom_indices.size();
		return;
	}

	void ComputeRemainingAttributes(){
		
		this->ComputeCentroid();	
		this->ComputeNormal();
		this->ComputeEffectiveRadius();
		return;
	}
	
	void ComputeCentroidAndNormal(){
		this->ComputeCentroid();
		this->ComputeNormal();
		return;
	}

	szv atom_indices;
	aptrv atom_ptrs;
	vptrv atom_coord_ptrs; 
	aptrv aromatic_atom_ptrs;
	vptrv aromatic_coord_ptrs;
	//For receptor atoms, these will point to atom.coords. For ligand atoms, this->coords[index].
	//Do this to avoid repetitively evaluating "if receptor atom, go here, if ligand atoms, go there". This slows down the scoring function.
	//atom.coords and model.coords, put receptor and ligand atom coordinate separately. This is a very stupid design. 
	//Or it could be smart, as an anti-repair mechanism done to the perfection. I mean, as smart as Apple.
	vec centroid;
	vec normal;
	fl effective_radius; //Average atom-centroid distance. 
};
typedef std::vector<ring_attribute> ring_info;

struct aliphatic_carbon_attribute {
	aliphatic_carbon_attribute(atom_vc* carbon_, vec* c_coord_, aptrv& h_neighbors_, vptrv& h_coords_, bool chpi_explicit_hydrogen){
		VINA_CHECK(h_neighbors_.size() == h_coords_.size());

		this->carbon = carbon_;
		this->c_coord = c_coord_;
		this->h_neighbors = h_neighbors_;
		this->h_coords = h_coords_;
		this->explicit_hydrogen = chpi_explicit_hydrogen;
	}

	void ComputeCHBondVectors(sz index, vec& bond){
		VINA_CHECK(this->explicit_hydrogen);
		VINA_CHECK(index < this->h_coords.size());

		vec* h_coord = this->h_coords[index];
		bond.data[0] = h_coord->data[0] - this->c_coord->data[0];
		bond.data[1] = h_coord->data[1] - this->c_coord->data[1];
		bond.data[2] = h_coord->data[2] - this->c_coord->data[2];
		return;
	}

	atom_vc* carbon;
	vec* c_coord;
	aptrv h_neighbors;
	vptrv h_coords;
	bool explicit_hydrogen;
};
typedef std::vector<aliphatic_carbon_attribute> aliphatic_carbon_info;

enum distance_type {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE};
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct non_cache; // forward declaration
struct naive_non_cache; // forward declaration
struct cache; // forward declaration
struct szv_grid; // forward declaration
struct terms; // forward declaration
struct conf_independent_inputs; // forward declaration
struct pdbqt_initializer; // forward declaration - only declared in parse_pdbqt.cpp
struct model_test;

struct model {
	void append(const model& m);
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }

	sz num_movable_atoms() const { return m_num_movable_atoms; }
	sz num_internal_pairs() const;
	sz num_other_pairs() const { return other_pairs.size(); }
	sz num_ligands() const { return ligands.size(); }
	sz num_flex() const { return flex.size(); }
	sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; }
	sz ligand_longest_branch(sz ligand_number) const;
	sz ligand_length(sz ligand_number) const;

	szv get_movable_atom_types(atom_type::t atom_typing_used_) const;

	conf_size get_size() const;
	conf get_initial_conf() const; // torsions = 0, orientations = identity, ligand positions = current

	grid_dims movable_atoms_box(fl add_to_each_dimension, fl granularity = 0.375) const;

	void write_flex  (                  const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }
	void write_ligand(sz ligand_number, const path& name, const std::string& remark) const { VINA_CHECK(ligand_number < ligands.size()); write_context(ligands[ligand_number].cont, name, remark); }
	void write_structure(ofile& out) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out);
		if(num_flex() > 0) // otherwise remark is written in vain
			write_context(flex_context, out);
	}
	void write_structure(ofile& out, const std::string& remark) const {
		out << remark;
		write_structure(out);
	}
	void write_structure(const path& name) const { ofile out(name); write_structure(out); }
	void write_model(ofile& out, sz model_number, const std::string& remark) const {
		out << "MODEL " << model_number << '\n';
		write_structure(out, remark);
		out << "ENDMDL\n";
	}
	void seti(const conf& c);
	void sete(const conf& c);
	void set (const conf& c);

	fl gyration_radius(sz ligand_number) const; // uses coords

	const atom_base& movable_atom  (sz i) const { assert(i < m_num_movable_atoms); return  atoms[i]; }
	const vec&       movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; }

	const vec& atom_coords(const atom_index& i) const;
	fl distance_sqr_between(const atom_index& a, const atom_index& b) const;
	bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const; // there is an atom closer to both a and b then they are to each other and immobile relative to them

	distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

	// clean up
	fl evali     (const precalculate& p,                  const vec& v                          ) const;
	fl evale     (const precalculate& p, const igrid& ig, const vec& v                          ) const;
	fl eval      (const precalculate& p, const igrid& ig, const vec& v, const conf& c           );
	double get_torsion_coords_vecs_list(vec A, vec B, vec C, vec D);
	double phi_alpha_energy(double phi_angle);
	double phi_beta_energy(double phi_angle);
	double psi_2A3E_energy(double psi_angle);
	double psi_2E3A_energy(double psi_angle);
	double psi_6A_energy(double psi_angle);
        double psi_6E_energy(double psi_angle);
        double omega_6A_energy(double omega_angle);
        double omega_6E_energy(double omega_angle);

	fl eval_deriv(const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g, const fl chi_coeff, const fl chi_cutoff);

	fl eval_intramolecular(                            const precalculate& p,                  const vec& v, const conf& c);
	fl eval_internal_torsional(const precalculate& p, const vec& v, const conf& c);
	fl eval_chi(const fl chi_coeff, const fl chi_cutoff);
	fl eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy);


	fl rmsd_lower_bound(const model& m) const; // uses coords
	fl rmsd_upper_bound(const model& m) const; // uses coords
	fl rmsd_ligands_upper_bound(const model& m) const; // uses coords

	void verify_bond_lengths() const;
	void about() const;

	vecv get_ligand_internal_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(internal_coords[i]);
		return tmp;
	}
	vecv get_flexible_coords();
	vecv get_ligand_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			{
			tmp.push_back(coords[i]);
			}
		return tmp;
	}
	vecv get_heavy_atom_movable_coords() const { // FIXME mv
		vecv tmp;
		VINA_FOR(i, num_movable_atoms())
			if(atoms[i].el != EL_TYPE_H)
				tmp.push_back(coords[i]);
		return tmp;
	}
	void check_internal_pairs() const;
	void print_stuff() const; // FIXME rm

	fl clash_penalty() const;
 
        bool bonded_to_HD(const atom_vc& a) const; //Yao: make this public, used to be private.
        bool bonded_to_heteroatom(const atom_vc& a) const; //Yao: make this public, used to be private. 
	bool bonded_to_hydrogen(const atom_vc& a) const;
	bool chpi_contact_possible(const atom_vc& central_carbon, vec& ring_centroid, vec& ring_normal, aptrv& ring, fl& ring_effective_radius, bool ring_is_receptor);

	void build_ar_ring_info(); //Yao added 2023602, find receptor aromatic rings, pre-compute ring centroids and normals. 
        void DFSVisit(aptrv& atoms, std::map<atom_vc*, bool>& atom_visited_map, std::map<atom_vc*, atom_vc*>& atom_parent_map, atom_vc* atom, int& counter, aptrv& current_path, std::vector<aptrv>& detected_cycles);

	fl eval_chpi();
	fl eval_chpi_c();
	fl eval_chpi_h();
	fl eval_chpi_c_ring(aliphatic_carbon_attribute& c_attr, ring_attribute& r_attr);
	fl eval_chpi_h_ring(aliphatic_carbon_attribute& c_attr, ring_attribute& r_attr);

	fl eval_chpi_enthalpy_c(fl& r);
	fl eval_chpi_enthalpy_h(fl& r);
	fl eval_chpi_entropy(fl& horizontal_offset);

	std::vector<aptrv> remove_duplicate_cycles(std::vector<aptrv>& cycles);
	std::vector<aptrv> remove_redundant_cycles(std::vector<aptrv>& cycles);

        std::vector<aptrv> DetectCyclesByDFS(aptrv& atoms);
        std::vector<aptrv> DetectAromaticCycles(aptrv& atoms);
	void DetectAromaticCycles(resv& residues, ring_info& ring_attributes);
	void build_residue_info();
	void build_residues_from_atoms(resv& residues, aptrv& atoms);

	//Yao added 20230602
        //Could contain multiple ligands. Each vector of aptrv is all the aromatic rings of that ligand.
	aliphatic_carbon_info lig_ali_carb_info;//cleared
	aliphatic_carbon_info rec_ali_carb_info;//cleared
	ring_info lig_ar_ring_info;//cleared
	ring_info rec_ar_ring_info;//cleared

	resv receptor_residues;//cleared
	resv ligand_residues;//cleared

	aptrv ligand_atoms;//cleared
	aptrv grid_atom_ptrs;//cleared
	fl weight_chpi;
	bool chpi_explicit_hydrogen;
	

private:
	friend struct non_cache;
	friend struct naive_non_cache;
	friend struct cache;
	friend struct szv_grid;
	friend struct terms;
	friend struct conf_independent_inputs;
	friend struct appender_info;
	friend struct pdbqt_initializer;
	friend struct model_test;
	friend struct parallel_mc_task; //Yao added 20230704

	model() : m_num_movable_atoms(0), m_atom_typing_used(atom_type::XS) {};

	const atom_vc& get_atom(const atom_index& i) const { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }
	      atom_vc& get_atom(const atom_index& i)       { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }

	void write_context(const context& c, ofile& out) const;
	void write_context(const context& c, ofile& out, const std::string& remark) const {
		out << remark;
	}
	void write_context(const context& c, const path& name) const {
		ofile out(name);
		write_context(c, out);
	}
	void write_context(const context& c, const path& name, const std::string& remark) const {
		ofile out(name);
		write_context(c, out, remark);
	}
	fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const; // actually static
	
	atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
	//bool bonded_to_HD(const atom_vc& a) const; //Yao: make it public for now
	//bool bonded_to_heteroatom(const atom_vc& a) const; //Yao: make it public for now
	sz find_ligand(sz a) const;
	void bonded_to(sz a, sz n, szv& out) const;
	szv bonded_to(sz a, sz n) const;

	void assign_bonds(const distance_type_matrix& mobility); // assign bonds based on relative mobility, distance and covalent length
	void assign_types();
	void initialize_pairs(const distance_type_matrix& mobility);
	void initialize(const distance_type_matrix& mobility);
	fl clash_penalty_aux(const interacting_pairs& pairs) const;

	vecv internal_coords;
	vecv coords;
	vecv minus_forces;

	atomv grid_atoms;
	atomv atoms; // movable, inflex
	vector_mutable<ligand> ligands;  //type declarations (?) 
	vector_mutable<residue_vc> flex;
	context flex_context;
	interacting_pairs other_pairs; // all except internal to one ligand: ligand-other ligands; ligand-flex/inflex; flex-flex/inflex

	sz m_num_movable_atoms;
	atom_type::t m_atom_typing_used;
};

#endif
