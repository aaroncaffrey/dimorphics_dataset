﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace dimorphics_dataset
{
    internal class atom
    {
        internal string uniprot_sequence;

        internal string pdb_id;
        internal int model_index;
        internal char chain_id;
        internal string atom_type;
        internal char element;
        internal char amino_acid;
        internal int residue_index = -1;
        internal int serial_index = -1;

        internal int array_index = -1;

        internal int master_index = -1;

        internal char i_code = ' ';

        internal char dssp_multimer = ' ';
        internal char dssp3_multimer = ' ';
        internal char dssp_monomer = ' ';
        internal char dssp3_monomer = ' ';

        internal char stride_multimer = ' ';
        internal char stride3_multimer = ' ';
        internal char stride_monomer = ' ';
        internal char stride3_monomer = ' ';

        // todo: unused variable: internal int is_strand_interface_atom;
        // todo: unused variable: internal int is_standard_strand_interface_atom;
        // todo: unused variable: internal int is_dimorphic_strand_interface_atom;
        // todo: unused variable: internal int is_hybrid_strand_interface_atom;

        internal char parallelism = '_';
        internal char symmetry_mode = '_';
        internal char strand_interface_type = '_';

        internal double X;
        internal double Y;
        internal double Z;

        internal info_solvent_access RSA_L;
        internal info_solvent_access RSA_S;

        internal info_sable_item sable_item;

        internal List<(string database, List<info_blast_pssm_entry> pssm_entries)> amino_acid_pssm_unnormalised;
        internal List<(string database, List<info_blast_pssm_entry> pssm_entries)> amino_acid_pssm_normalised;

        internal List<(string format, info_mpsa_reader_line_entry mpsa_entry)> mpsa_entries;

        internal (int index, double short_type_score, double long_type_score, double glob_type_score, double anchor2_score) iup_entry;

        internal List<atom> amino_acid_atoms;
        internal List<atom> chain_atoms;

        internal atom amino_acid_master_atom;
        //internal List<(atom atom, double distance)> contact_map_intermolecular;
        //internal List<(atom atom, double distance)> contact_map_intramolecular;

        internal List<(int array_index, char chain_id)> chain_order_in_pdb_file;

        internal double foldx_monomer_ala_scan_ddg;

        internal List<info_ring_edge> monomer_ring_edges;
        internal List<info_ring_node> monomer_ring_nodes;


        // todo: unused variable: internal bool amino_acid_exists_in_repaired_structure;

        internal double chain_dna_binding_prob_nr;
        internal double chain_dna_binding_prob_swissprot;
        internal double chain_dna_binding_prob_uniref90;

        //internal List<(bool repaired, char mutant_res, double ddg)> foldx_position_scan_ddg;

        //META30P

        internal (int index1, int index2, double distance)[] intramolecular_contact_flat_table;
        // todo: unused variable: internal (atom atom1, atom atom2, double distance)[] intramolecular_contact_flat_ref_table;
        internal double[][] intramolecular_contact_table;
        internal int intramolecular_contact_table_index = -1;

        internal static double[][] make_intramolecular_contact_table(List<atom> p)
        {
            if (p == null)
            {
                return null;
            }

            var d_flat = new (int index1, int index2, double distance)[p.Count * p.Count]; // could use only half ... ((p.Count*p.Count)/2)-(p.Count/2)
            var d_ref_flat = new (atom atom1, atom atom2, double distance)[p.Count * p.Count];
            var d = new double[p.Count][];

            for (var i = 0; i < p.Count; i++)
            {
                d[i] = new double[p.Count];
            }

            var k = -1;
            for (var i = 0; i < p.Count; i++)
            {
                for (var j = 0; j < p.Count; j++)
                {
                    k++;
                    if (i == j)
                    {
                        d_flat[k] = (i, j, 0d);
                        d_flat[k] = (j, i, 0d);
                        d_ref_flat[k] = (p[i], p[j], 0d);
                        d_ref_flat[k] = (p[j], p[i], 0d);
                    }

                    if (i < j) continue;

                    var x = atom.Distance3D(p[i], p[j]);
                    d[i][j] = x;
                    d[j][i] = x;
                    d_flat[k] = (i, j, x);
                    d_flat[k] = (j, i, x);
                    d_ref_flat[k] = (p[i], p[j], x);
                    d_ref_flat[k] = (p[j], p[i], x);
                }
            }

            for (var i = 0; i < p.Count; i++)
            {
                p[i].intramolecular_contact_table_index = i;
                p[i].intramolecular_contact_table = d;
                p[i].intramolecular_contact_flat_table = d_flat;
            }

            return d;
        }

        internal static char Aa3To1(string aa)
        {
            var aa_codes = new (string three_letter_code, char one_letter_code)[]
            {
                ("Arg", 'R'), ("Lys", 'K'), ("Asp", 'D'), ("Glu", 'E'), ("Gln", 'Q'), ("Asn", 'N'), ("His", 'H'),
                ("Ser", 'S'), ("Thr", 'T'), ("Tyr", 'Y'), ("Cys", 'C'), ("Met", 'M'), ("Trp", 'W'), ("Ala", 'A'),
                ("Ile", 'I'), ("Leu", 'L'), ("Phe", 'F'), ("Val", 'V'), ("Pro", 'P'), ("Gly", 'G'), ("XAA", 'X'),
                ("UNK", 'X'), ("MSE", 'M')
            };

            var r = aa_codes.FirstOrDefault(a => String.Equals(a.three_letter_code, aa.Trim(), StringComparison.OrdinalIgnoreCase)).one_letter_code;

            if (r == default(char))
            {
            }

            return r;
        }

        internal static string Aa1To3(char aa)
        {
            if (!char.IsLetter(aa)) return /*program.string_debug*/($@"{aa}");

            var aa_codes = new (string three_letter_code, char one_letter_code)[]
            {
                ("Arg", 'R'), ("Lys", 'K'), ("Asp", 'D'), ("Glu", 'E'), ("Gln", 'Q'), ("Asn", 'N'), ("His", 'H'),
                ("Ser", 'S'), ("Thr", 'T'), ("Tyr", 'Y'), ("Cys", 'C'), ("Met", 'M'), ("Trp", 'W'), ("Ala", 'A'),
                ("Ile", 'I'), ("Leu", 'L'), ("Phe", 'F'), ("Val", 'V'), ("Pro", 'P'), ("Gly", 'G'), ("XAA", 'X'),
                ("UNK", 'X')
            };

            var r = aa_codes.FirstOrDefault(a => a.one_letter_code == aa).three_letter_code;

            if (r == default)
            {
                r = /*program.string_debug*/($@"{aa}");
            }

            return r.ToUpperInvariant();
        }

        internal static string secondary_structure_state_reduction(string ss)
        {
            if (string.IsNullOrEmpty(ss)) return ss;

            return string.Concat(ss.Select(secondary_structure_state_reduction).ToList());
        }

        internal static char secondary_structure_state_reduction(char ss)
        {
            if (ss == 'E' || ss == 'B' || ss == 'e' || ss == 'b') return 'E'; // change E/e, B/b to E
            if (ss == 'G' || ss == 'g' || ss == 'H' || ss == 'h') return 'H'; // change G/g, H/h to H
            //if (ss == 'T' || ss == 'S' || ss == ' ' || ss == '_') return 'C';
            return 'C'; // all others to C
        }

        internal static string[] extract_split_pdb_chains(string pdb_id, char? chain_id)
        {
            var use_ram_disk = false;

            var drive_letter = use_ram_disk ? /*program.string_debug*/($@"r") : /*program.string_debug*/($@"c");

            var pdb_in_folder = /*program.string_debug*/($@"{drive_letter}:\phd\betastrands_dataset\pdb\");
            var pdb_out_folder = /*program.string_debug*/($@"{drive_letter}:\phd\betastrands_dataset\pdb_split\");

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            pdb_id = pdb_id.ToUpperInvariant();

            var pdb_lines = io_proxy.ReadAllLines(/*program.string_debug*/($@"{pdb_in_folder}{pdb_id}.pdb"), nameof(atom), nameof(extract_split_pdb_chains)).ToList();
            var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith(/*program.string_debug*/($@"ENDMDL"), StringComparison.Ordinal));
            if (endmdl_index > -1) pdb_lines = pdb_lines.Take(endmdl_index).ToList();
            pdb_lines = pdb_lines.Where(a => a.StartsWith(/*program.string_debug*/($@"ATOM "), StringComparison.Ordinal)).ToList();

            var chains = pdb_lines.GroupBy(a => a[21]).Where(a => chain_id == null || a.Key == chain_id).ToList();
            chains.ForEach(a =>
            {
                var chain_pdb_file = Path.Combine(pdb_out_folder, /*program.string_debug*/($@"{pdb_id}{a.Key}.pdb"));

                if (!File.Exists(chain_pdb_file) || new FileInfo(chain_pdb_file).Length == 0)
                {
                    //Directory.CreateDirectory(Path.GetDirectoryName(chain_pdb_file));

                    io_proxy.WriteAllLines(chain_pdb_file, a.ToList(), nameof(atom), nameof(extract_split_pdb_chains));
                }
            });
            return chains.Select(a => /*program.string_debug*/($@"{pdb_out_folder}{pdb_id}{a.Key}.pdb")).ToArray();
        }

        internal static string extract_split_pdb_chains_res_ids(string pdb_id, char chain_id, int first_res_id, int last_res_id)
        {
            if (string.IsNullOrWhiteSpace(pdb_id))
            {
                throw new ArgumentNullException(nameof(pdb_id));
            }

            if (first_res_id > last_res_id)
            {
                var first = last_res_id;
                var last = first_res_id;
                first_res_id = first;
                last_res_id = last;
            }

            var use_ram_disk = false;

            var drive_letter = use_ram_disk ? /*program.string_debug*/($@"r") : /*program.string_debug*/($@"c");

            var pdb_in_folder = /*program.string_debug*/($@"{drive_letter}:\phd\betastrands_dataset\pdb\");
            var pdb_out_folder = /*program.string_debug*/($@"{drive_letter}:\phd\betastrands_dataset\pdb_split\");

            pdb_id = pdb_id.ToUpperInvariant();
            var pdb_lines = io_proxy.ReadAllLines(/*program.string_debug*/($@"{pdb_in_folder}{pdb_id}.pdb"), nameof(atom), nameof(extract_split_pdb_chains_res_ids)).ToList();
            var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith(/*program.string_debug*/($@"ENDMDL"), StringComparison.Ordinal));
            if (endmdl_index > -1) pdb_lines = pdb_lines.Take(endmdl_index).ToList();
            pdb_lines = pdb_lines.Where(a => a.StartsWith(/*program.string_debug*/($@"ATOM "), StringComparison.Ordinal) && a[21] == chain_id && int.Parse(a.Substring(22, 4), NumberStyles.Integer, NumberFormatInfo.InvariantInfo) >= first_res_id && int.Parse(a.Substring(22, 4), NumberStyles.Integer, NumberFormatInfo.InvariantInfo) <= last_res_id).ToList();

            //this.residue_index = int.Parse(pdb_atom_line.Substring(22, 4));
            //this.i_code = pdb_atom_line[26];

            var output_pdb_file = Path.Combine(/*program.string_debug*/($@"{pdb_out_folder}"), /*program.string_debug*/($@"{pdb_id}{chain_id}_{first_res_id}_{last_res_id}.pdb"));

            if (!File.Exists(output_pdb_file) || new FileInfo(output_pdb_file).Length == 0)
            {
                //Directory.CreateDirectory(Path.GetDirectoryName(output_pdb_file));

                io_proxy.WriteAllLines(output_pdb_file, pdb_lines.ToList(), nameof(atom), nameof(extract_split_pdb_chains_res_ids));
            }

            return output_pdb_file;
        }

        internal static void run_dssp(string pdb_id, string pdb_folder = @"pdb", string dssp_exe = @"dssp2.exe")
        {
            pdb_folder = /*program.string_debug*/(Path.Combine(program.data_root_folder, pdb_folder));
            dssp_exe = /*program.string_debug*/(Path.Combine(program.data_root_folder, pdb_folder, dssp_exe));

            var pdb_file = /*program.string_debug*/($@"{pdb_folder}{pdb_id}.pdb");
            var dssp_file = /*program.string_debug*/($@"{pdb_folder}{pdb_id}.dssp");


            var start = new ProcessStartInfo
            {
                WorkingDirectory = Path.GetDirectoryName(dssp_exe) ?? /*program.string_debug*/($@""),
                FileName = dssp_exe,
                Arguments = /*program.string_debug*/($@"-i ""{pdb_file}"" -o ""{dssp_file}"""),
                UseShellExecute = false,
                CreateNoWindow = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true
            };

            using var process = Process.Start(start);
            if (process == null) return;

            io_proxy.WriteLine(/*program.string_debug*/($@"{start.FileName} {start.Arguments}"), nameof(atom), nameof(run_dssp));
            using var reader = process.StandardOutput;
            var result = reader.ReadToEnd();
            if (!string.IsNullOrWhiteSpace(result)) io_proxy.WriteLine(result, nameof(atom), nameof(run_dssp));

            var stderr = process.StandardError.ReadToEnd();
            if (!string.IsNullOrWhiteSpace(stderr)) io_proxy.WriteLine(stderr, nameof(atom), nameof(run_dssp));
            return;
        }

        internal atom(string pdb_id, string pdb_atom_line, int pdb_model_array_index = -1, int array_index = -1)
        {
            if (string.IsNullOrWhiteSpace(pdb_atom_line))
            {
                throw new ArgumentNullException(pdb_atom_line);
            }

            this.pdb_id = pdb_id;
            this.chain_id = pdb_atom_line[21];
            this.residue_index = int.Parse(pdb_atom_line.Substring(22, 4), NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            this.i_code = pdb_atom_line[26];
            this.amino_acid = atom.Aa3To1(pdb_atom_line.Substring(17, 3));
            this.array_index = array_index;
            this.model_index = pdb_model_array_index;

            this.serial_index = int.Parse(pdb_atom_line.Substring(7 - 1, (11 - 7) + 1), NumberStyles.Integer, NumberFormatInfo.InvariantInfo);

            this.atom_type = pdb_atom_line.Substring(13 - 1, (16 - 13) + 1).Trim(); // is this CA / CB  etc ??
            this.element = pdb_atom_line.Length >= 77 - 1 ? pdb_atom_line.Substring(77 - 1, (78 - 77) + 1).Trim().FirstOrDefault() : default;
            this.X = double.Parse(pdb_atom_line.Substring(31 - 1, (38 - 31) + 1), NumberStyles.Float, NumberFormatInfo.InvariantInfo);
            this.Y = double.Parse(pdb_atom_line.Substring(39 - 1, (46 - 39) + 1), NumberStyles.Float, NumberFormatInfo.InvariantInfo);
            this.Z = double.Parse(pdb_atom_line.Substring(47 - 1, (54 - 47) + 1), NumberStyles.Float, NumberFormatInfo.InvariantInfo);
        }

        internal atom(string pdb_id, char chain_id, int residue_index, char i_code, char amino_acid, int array_index, int model_index, int serial_index, string atom_type, char element, double X, double Y, double Z)
        {
            this.pdb_id = pdb_id;
            this.chain_id = chain_id;
            this.residue_index = residue_index;
            this.i_code = i_code;
            this.amino_acid = amino_acid;
            this.array_index = array_index;
            this.model_index = model_index;
            this.serial_index = serial_index;
            this.atom_type = atom_type;
            this.element = element;
            this.X = X;
            this.Y = Y;
            this.Z = Z;
        }

        //internal static bool use_ram_disk = true;
        //internal const string pdb_folder = /*/*program.string_debug*/($@"{(use_ram_disk?"r":"c")}*/ Path.Combine(program.data_root_folder,"pdb\";

        internal static List<(string pdb_id, int pdb_model_index, char chain_id, List<atom> pdb_model_chain_atoms)> load_atoms_pdb
        (
            string pdb_id,

            load_atoms_pdb_options options,

            string pdb_folder = null
        )
        {
            if (options == null)
            {
                throw new ArgumentNullException(nameof(options));
            }

            if (string.IsNullOrWhiteSpace(pdb_folder))
            {
                pdb_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"pdb"));
            }


            //var load_ring_data = true;
            //var load_mpsa_sec_struct_predictions = true;
            //var load_blast_pssms = true;
            //var load_iup_data = true;
            //var load_sable = true;
            //var load_ala_scan = true;
            //var load_dna_binding_vars = true;

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);
            pdb_id = pdb_id.ToUpperInvariant();
            var pdb_id_simple = pdb_id.Substring(0, 4);
            var result = new List<(string pdb_id, int pdb_model_index, char chain_id, List<atom> pdb_model_chain_atoms)>();
            var pdb_filename = Path.Combine(pdb_folder, /*program.string_debug*/($@"{pdb_id}.pdb"));
            var pdb_lines = io_proxy.ReadAllLines(pdb_filename, nameof(atom), nameof(load_atoms_pdb)).ToList();
            var pdb_model_array_index = -1;
            var file_cache = new List<(string filename, string[] data)>();

            while (pdb_lines.Count > 0)
            {
                // set model index
                pdb_model_array_index++;

                // find end of model index
                var endmdl_index = pdb_lines.FindIndex(a => a.StartsWith(/*program.string_debug*/($@"ENDMDL"), StringComparison.Ordinal));

                // get all pdb lines for this model
                var pdb_model_lines = (endmdl_index > -1) ? pdb_lines.Take(endmdl_index + 1).ToList() : pdb_lines;

                // remove model from pdb_lines
                pdb_lines = pdb_lines.Skip(pdb_model_lines.Count).ToList();

                // filter model for ATOMS only
                pdb_model_lines = pdb_model_lines.Where(a => a.StartsWith(/*program.string_debug*/($@"ATOM "), StringComparison.Ordinal)).ToList();

                // if no ATOMs, go to next model
                if (pdb_model_lines.Count == 0) continue;

                // convert ATOM records to ATOM classes
                var pdb_model_atoms = pdb_model_lines.Select((atom_struct, array_index) => new atom(pdb_id, atom_struct, pdb_model_array_index, array_index)).ToList();

                // save the chain order, which may not be alphabetical, or such
                var chain_order_in_pdb_file = pdb_model_atoms.Select(atom => atom.chain_id).Distinct().Select((chain_id, array_index) => (array_index, chain_id)).ToList();
                pdb_model_atoms.ForEach(a => a.chain_order_in_pdb_file = chain_order_in_pdb_file);

                // re-order the atom classes by chain id, res_id, i_code
                pdb_model_atoms = pdb_model_atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

                // remove any additional i_codes
                if (options.first_icode_only)
                {
                    pdb_model_atoms = pdb_model_atoms.GroupBy(a => (a.chain_id, a.residue_index)).SelectMany(a => a.ToList().GroupBy(b => b.i_code).First().ToList()).ToList();
                }

                // make list of atoms for each amino acid
                foreach (var atom in pdb_model_atoms)
                {
                    if (atom.amino_acid_atoms != null) continue;
                    atom.amino_acid_atoms = pdb_model_atoms.Where(a => a.model_index == atom.model_index && string.Equals(a.pdb_id, atom.pdb_id, StringComparison.OrdinalIgnoreCase) && a.chain_id == atom.chain_id && a.residue_index == atom.residue_index && a.i_code == atom.i_code).ToList();
                    atom.amino_acid_atoms.ForEach(a => a.amino_acid_atoms = atom.amino_acid_atoms);
                }

                // get the main/master atoms
                var pdb_model_master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);

                // set the master_atom for all atoms
                for (var master_atom_index = 0; master_atom_index < pdb_model_master_atoms.Count; master_atom_index++)
                {
                    var master_atom = pdb_model_master_atoms[master_atom_index];

                    master_atom.amino_acid_atoms.ForEach(a => a.amino_acid_master_atom = master_atom);
                }

                // load ring values
                if (options.load_3d_ring_data) { atom.load_ring(pdb_id, pdb_model_atoms); }

                // load dssp values
                if (options.load_3d_dssp_data) { atom.load_dssp(pdb_id, pdb_model_atoms); }

                // load stride values
                if (options.load_3d_stride_data) { atom.load_stride(pdb_id, pdb_model_atoms); }

                // load free-sasa values (rsa)
                if (options.load_3d_rsa_data) { atom.load_rsa(pdb_id, pdb_model_atoms); }


                // load psi-blast PSSMs

                if (options.load_1d_blast_pssms) atom.load_pssm(pdb_id, pdb_model_atoms);

                // load MPSA secondary structure predictions
                if (options.load_2d_mpsa_sec_struct_predictions)
                {
                    atom.load_mpsa_sec_struct_predictions(pdb_id, pdb_model_atoms);

                    //pdb_model_chain_atoms.ForEach(a => compare_dssp_to_mpsa(pdb_id, a.chain_id, a.pdb_model_chain_atoms));
                }


                // find intra/inter-molecular contacts/interactions


                //if (options.find_3d_intermolecular)
                //{
                //    atom.find_intermolecular_contacts(pdb_id, pdb_model_atoms, (double) 5);
                //}

                //if (options.find_3d_intramolecular)
                //{
                //    atom.find_intramolecular_contacts(pdb_id, pdb_model_atoms, (double) 5);
                //}


                // further split model atoms into chain atoms
                var pdb_model_chain_atoms = pdb_model_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, pdb_model_chain_atoms: a.ToList())).ToList();

                if (options.find_3d_intramolecular)
                {
                    foreach (var c in pdb_model_chain_atoms)
                    {
                        var d = make_intramolecular_contact_table(c.pdb_model_chain_atoms);
                    }
                }

                var pdb_model_chain_master_atoms = pdb_model_chain_atoms.Select(a => (chain_id: a.chain_id, pdb_model_chain_master_atoms: select_amino_acid_master_atoms(pdb_id, a.pdb_model_chain_atoms))).ToList();

                foreach (var c in pdb_model_chain_master_atoms)
                {
                    // get the main/master atoms
                    //var pdb_model_chain_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, c.pdb_model_chain_atoms);


                    // set the master_atom for all atoms
                    for (var master_atom_index = 0; master_atom_index < c.pdb_model_chain_master_atoms.Count; master_atom_index++)
                    {
                        var master_atom = c.pdb_model_chain_master_atoms[master_atom_index];

                        master_atom.amino_acid_atoms.ForEach(a => a.amino_acid_master_atom = master_atom);
                        master_atom.amino_acid_atoms.ForEach(a => a.master_index = master_atom_index);
                    }
                }

                if (options.load_1d_iup_data)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        // get the main/ master atoms
                        //var pdb_model_chain_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, c.pdb_model_chain_atoms);

                        //var chain_seq = string.Join(/*program.string_debug*/($@""), pdb_model_chain_master_atoms.Select(a => a.amino_acid).ToList());


                        var iup_data = info_iupred.load(pdb_id_simple, c.chain_id);

                        if (iup_data == null || iup_data.Count == 0)
                        {
                            continue;
                        }

                        for (var master_atom_index = 0; master_atom_index < c.pdb_model_chain_master_atoms.Count; master_atom_index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[master_atom_index];

                            master_atom.iup_entry = iup_data[master_atom_index];

                            foreach (var atom in c.pdb_model_chain_master_atoms[master_atom_index].amino_acid_atoms)
                            {
                                atom.iup_entry = master_atom.iup_entry;
                            }
                        }
                    }
                }

                if (options.load_1d_dna_binding)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        load_dna_binding_stackdppred(pdb_id + c.chain_id, c.pdb_model_chain_master_atoms);
                    }
                }

                if (options.load_3d_foldx_ala_scan)
                {
                    foreach (var c in pdb_model_chain_atoms)
                    {
                        // load foldx AlaScan values

                        // 1EV0A_Repair_AS.fxout
                        var alascan = info_foldx.load_foldx_ala_scanning(/*program.string_debug*/($@"{pdb_id}{c.chain_id}_Repair"), c.chain_id, null, false).data;

                        if (alascan != null && alascan.Count > 0)
                        {
                            foreach (var scan in alascan)
                            {
                                var atoms = c.pdb_model_chain_atoms.Where(a => string.Equals(/*program.string_debug*/($@"{a.pdb_id}{a.chain_id}_Repair"), scan.pdb_id, StringComparison.OrdinalIgnoreCase) && a.chain_id == scan.chain_id && a.residue_index == scan.residue_index).ToList();
                                atoms.ForEach(a => a.foldx_monomer_ala_scan_ddg = scan.ddg);
                            }
                        }
                    }

                    // load foldx PositionScan values
                    // var load_pos_scan = true;
                    // if (load_pos_scan) Atom.load_foldx_position_scanning(pdb_model_atoms);
                }

                if (options.load_uniprot)
                {
                    // load uniprot sequence (& still todo: also load mapping and save in atoms)
                    //var chain_ids = pdb_model_chain_atoms.Select(a => a.chain_id).Distinct().ToList();
                    foreach (var c in pdb_model_chain_atoms)
                    //    foreach (var chain_id in chain_ids)
                    {
                        var uniprot_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"uniprot"),
                            /*program.string_debug*/($@"{pdb_id_simple}{c.chain_id}.fasta"));

                        var uniprot_file_data = file_cache.FirstOrDefault(a => string.Equals(a.filename, uniprot_file, StringComparison.Ordinal)).data;

                        if (uniprot_file_data == null)
                        {
                            if (!File.Exists(uniprot_file) || new FileInfo(uniprot_file).Length == 0) continue;
                            uniprot_file_data =
                                io_proxy.ReadAllLines(uniprot_file, nameof(atom), nameof(load_atoms_pdb));
                            file_cache.Add((uniprot_file, uniprot_file_data));
                        }

                        var uniprot_sequence = string.Join(/*program.string_debug*/($@""),
                            uniprot_file_data.Where(a =>
                                    !string.IsNullOrWhiteSpace(a) &&
                                    !a.StartsWith(/*program.string_debug*/($@">"), StringComparison.Ordinal))
                                .ToList());

                        c.pdb_model_chain_atoms.ForEach(a => a.uniprot_sequence = uniprot_sequence);
                    }
                }

                // load sable
                if (options.load_1d_sable)
                {
                    foreach (var c in pdb_model_chain_master_atoms)
                    {
                        var sable_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"sable"), /*program.string_debug*/($@"{pdb_id_simple}{c.chain_id}.txt"));
                        var sable_data = info_sable.load(sable_file);


                        // check sable sequences matches this pdb
                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            c.pdb_model_chain_master_atoms[index].sable_item = sable_data[index];

                            if (c.pdb_model_chain_master_atoms[index].amino_acid != sable_data[index].amino_acid) throw new Exception();
                        }

                        // load sable data into atoms
                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[index];

                            master_atom.sable_item = sable_data[index];
                        }

                        // load ss as mpsa instance
                        var list = sable_data.Select(a => (a.seq_index, a.amino_acid, a.predicted_ss, a.prob_h, a.prob_e, a.prob_e)).ToList();
                        var mpsa_reader_sable = new info_mpsa_reader(/*program.string_debug*/($@"sable"), list);

                        for (var index = 0; index < c.pdb_model_chain_master_atoms.Count; index++)
                        {
                            var master_atom = c.pdb_model_chain_master_atoms[index];
                            if (master_atom.mpsa_entries == null) master_atom.mpsa_entries = new List<(string format, info_mpsa_reader_line_entry)>();
                            var mpsa_matrix_line_entry = mpsa_reader_sable?.mpsa_matrix?.First /*OrDefault*/(a => a.index == index);


                            master_atom.mpsa_entries.Add((mpsa_reader_sable.format, mpsa_matrix_line_entry));

                            foreach (var atom in master_atom.amino_acid_atoms)
                            {
                                atom.mpsa_entries = master_atom.mpsa_entries;
                            }
                        }
                    }
                }

                // add split model/chain to returned result
                pdb_model_chain_atoms.ForEach(a => result.Add((pdb_id, pdb_model_array_index, a.chain_id, a.pdb_model_chain_atoms)));


                // set chain_atoms property
                pdb_model_chain_atoms.ForEach(a => a.pdb_model_chain_atoms.ForEach(b => b.chain_atoms = a.pdb_model_chain_atoms));

                // output data for comparison of DSSP & STRIDE data to MPSA secondary structure predictions
                //pdb_model_chain_atoms.ForEach(a => compare_dssp_to_mpsa(pdb_id, a.chain_id, a.pdb_model_chain_atoms));

                if (options.first_model_only)
                {
                    break;
                }
            }

            if (options.first_model_only)
            {
                result = result.Where(a => a.pdb_model_index == 0).ToList();
            }

            if (options.load_3d_dssp_data && options.load_2d_mpsa_sec_struct_predictions)
            {
                foreach (var r in result)
                {
                    var dssp_mpsa = get_dssp_and_mpsa_subsequences(/*r.pdb_id, r.chain_id, */r.pdb_model_chain_atoms, (enum_get_dssp_and_mpsa_subsequences_params)0b_1111_1111_1111);

                    var ground_truths = dssp_mpsa; // dssp_mpsa.Where(a => a.format.StartsWith(/*program.string_debug*/($@"monomer", StringComparison.OrdinalIgnoreCase) || a.format.StartsWith(/*program.string_debug*/($@"multimer", StringComparison.OrdinalIgnoreCase)).ToList();

                    var data = new List<string>();

                    var header = /*program.string_debug*/($@"pdb,chain_id,format,prediction,{string.Join(/*program.string_debug*/($@","), ground_truths.Select(a => a.format).ToList())}");

                    data.Add(header);

                    foreach (var a in dssp_mpsa)
                    {
                        var x = /*program.string_debug*/($@"{r.pdb_id},{r.chain_id},{a.format},{a.prediction}");

                        foreach (var g in ground_truths)
                        {
                            var score = 0;
                            for (var i = 0; i < g.prediction.Length; i++)
                            {
                                if (g.prediction[i] == a.prediction[i]) score++;
                            }

                            var p = g.prediction.Length != 0 ? (double)score / (double)g.prediction.Length : 0d;

                            x = /*program.string_debug*/($@"{x},{p}");
                        }

                        data.Add(x);
                    }

                    data.Add(/*program.string_debug*/($@""));

                    lock (dssp_mpsa_lock)
                    {
                        var q3_fn = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"dssp_mpsa_protein_q3_data.csv"));
                        //Directory.CreateDirectory(Path.GetDirectoryName(q3_fn));
                        io_proxy.AppendAllLines(q3_fn, data, nameof(atom), nameof(load_atoms_pdb));
                    }
                }
            }

            return result;
        }

        private static readonly object dssp_mpsa_lock = new object();

        private static readonly List<string> compare_dssp_to_mpsa_list = new List<string>();

        private static readonly object compare_dssp_to_mpsa_lock = new object();

        internal static void compare_dssp_to_mpsa(string pdb_id, char chain_id, List<atom> param_atoms)
        {
            lock (compare_dssp_to_mpsa_lock)
            {
                if (compare_dssp_to_mpsa_list.Contains(/*program.string_debug*/($@"{pdb_id}{chain_id}"))) return;
                compare_dssp_to_mpsa_list.Add(/*program.string_debug*/($@"{pdb_id}{chain_id}"));

                var result1 = new List<(string filename, string pdb, string format, string data_category1, string data_category2, string[] data)>();
                var result2 = new List<(string filename, string pdb, string format, string data_category1, string data_category2, string[] data)>();

                var pdb = /*program.string_debug*/($@"{pdb_id}{chain_id}");

                var master_atoms = atom.select_amino_acid_master_atoms(null, param_atoms);

                var filename = /*program.string_debug*/($@"");
                var format = /*program.string_debug*/($@"pdb");

                var aa_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.amino_acid).ToArray());
                var aa_seq_aa_types = aa_seq.Distinct().OrderBy(a => a).ToArray();

                var strand_interface_type_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.strand_interface_type).ToArray());
                var strand_interface_type_seq_aa_types = strand_interface_type_seq.Distinct().OrderBy(a => a).ToArray();

                var parallelism_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.parallelism).ToArray());
                var parallelism_seq_aa_types = parallelism_seq.Distinct().OrderBy(a => a).ToArray();

                var symmetry_mode_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.symmetry_mode).ToArray());
                var symmetry_mode_seq_aa_types = symmetry_mode_seq.Distinct().OrderBy(a => a).ToArray();

                var dssp_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp_monomer).ToArray());
                var dssp_seq_ss_types = dssp_seq.Distinct().OrderBy(a => a).ToArray();

                var stride_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride_monomer).ToArray());
                var stride_seq_ss_types = stride_seq.Distinct().OrderBy(a => a).ToArray();

                var dssp3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp3_monomer).ToArray());
                var dssp3_seq_ss_types = dssp3_seq.Distinct().OrderBy(a => a).ToArray();

                var stride3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride3_monomer).ToArray());
                var stride3_seq_ss_types = stride3_seq.Distinct().OrderBy(a => a).ToArray();

                
                // add sequences to result
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_index"), /*program.string_debug*/($@"seq_index"), aa_seq.Select((b, i) => (i + 1).ToString(CultureInfo.InvariantCulture)).ToArray()));

                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_aa"), /*program.string_debug*/($@"seq_aa"), aa_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_bsi"), /*program.string_debug*/($@"seq_bsi"), strand_interface_type_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_symmetry"), /*program.string_debug*/($@"seq_symmetry"), symmetry_mode_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_parallel"), /*program.string_debug*/($@"seq_parallel"), parallelism_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_dssp"), /*program.string_debug*/($@"seq_dssp"), dssp_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_dssp3"), /*program.string_debug*/($@"seq_dssp3"), dssp3_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_stride"), /*program.string_debug*/($@"seq_stride"), stride_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_stride3"), /*program.string_debug*/($@"seq_stride3"), stride3_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));

                // add numeric values to result
                //aa_seq_aa_types.ToList().ForEach(aa => result.Add((filename, pdb, format, nameof(aa_seq), aa.ToString(CultureInfo.InvariantCulture), aa_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //strand_interface_type_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, /*program.string_debug*/($@"seq_bsi", aa.ToString(CultureInfo.InvariantCulture), strand_interface_type_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //symmetry_mode_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, /*program.string_debug*/($@"seq_symmetry", aa.ToString(CultureInfo.InvariantCulture), symmetry_mode_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));                
                //parallelism_seq_aa_types.Where(a => a != '_').ToList().ForEach(aa => result2.Add((filename, pdb, format, /*program.string_debug*/($@"seq_parallel", aa.ToString(CultureInfo.InvariantCulture), parallelism_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));                
                //dssp_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(dssp_seq), aa.ToString(CultureInfo.InvariantCulture), dssp_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //dssp3_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(dssp3_seq), aa.ToString(CultureInfo.InvariantCulture), dssp3_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //stride_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(stride_seq), aa.ToString(CultureInfo.InvariantCulture), stride_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //stride3_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, nameof(stride3_seq), aa.ToString(CultureInfo.InvariantCulture), stride3_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));
                //var atoms_all_dssp_seq = string.Join(/*program.string_debug*/($@""), atoms_all.Select(a => a.monomer_dssp).ToArray());
                //var atoms_all_dssp_seq_ss_types = atoms_all_dssp_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_stride_seq = string.Join(/*program.string_debug*/($@""), atoms_all.Select(a => a.stride_monomer).ToArray());
                //var atoms_all_stride_seq_ss_types = atoms_all_stride_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_dssp3_seq = string.Join(/*program.string_debug*/($@""), atoms_all.Select(a => a.monomer_dssp3).ToArray());
                //var atoms_all_dssp3_seq_ss_types = atoms_all_dssp3_seq.Distinct().OrderBy(a => a).ToArray();
                //var atoms_all_stride3_seq = string.Join(/*program.string_debug*/($@""), atoms_all.Select(a => a.monomer_stride3).ToArray());
                //var atoms_all_stride3_seq_ss_types = atoms_all_stride3_seq.Distinct().OrderBy(a => a).ToArray();

                var mpsa_seqs = master_atoms.SelectMany(a => a.mpsa_entries).ToArray();

                var q3 = new List<(string pdb_id, char chain_id, string truth_format, string predictor_format, string aa_subset, char ss, double q3_value, double truth_total)> { (/*program.string_debug*/($@"pdb_id"), 'c', /*program.string_debug*/($@"truth_format"), /*program.string_debug*/($@"predictor_format"), /*program.string_debug*/($@"aa_subset"), 's', 0, 0) };

                // for each subset, what is the average predicted value
                var av = new List<(string pdb_id, char chain_id, string predictor_format, string aa_subset, char ss, double average, double truth_total)> { (/*program.string_debug*/($@"pdb_id"), 'c', /*program.string_debug*/($@"predictor_format"), /*program.string_debug*/($@"aa_subset"), 's', 0, 0) };

                // problem(?): this doesn't separate the interfaces

                //const string  = ;

                foreach (var s in mpsa_seqs.GroupBy(a => a.format))
                {
                    format = s.Key;

                    {
                        // all
                        var atoms_all = master_atoms.Where(a => true).ToList();
                        var atoms_all_mpsa = atoms_all.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_all_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_all.Select((a, i) => ((m == '*' || atoms_all[i].dssp3_monomer == m || atoms_all[i].dssp_monomer == m) && (atoms_all[i].dssp3_monomer == atoms_all_mpsa[i].predicted_ss_code || atoms_all[i].dssp_monomer == atoms_all_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_all.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_all_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_all), a.ss, a.value, a.truth_total)));
                        var atoms_all_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_all.Select((a, i) => ((m == '*' || atoms_all[i].stride3_monomer == m || atoms_all[i].stride_monomer == m) && (atoms_all[i].stride3_monomer == atoms_all_mpsa[i].predicted_ss_code || atoms_all[i].stride_monomer == atoms_all_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_all.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_all_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_all), a.ss, a.value, a.truth_total)));
                        var atoms_all_mpsa_average = atoms_all_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList(); // each SS type, average probability
                        atoms_all_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_all), a.ss, a.av, a.truth_total)));
                        /*program.string_debug*/($@"HEC").ToList().ForEach(a => av.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), nameof(atoms_all), a, atoms_all.Count > 0 ? (double)atoms_all.Count(b => b.dssp3_monomer == a || b.dssp_monomer == a) / (double)atoms_all.Count : 0d, atoms_all.Count))); //"HECT"
                    }

                    {
                        // non-interface
                        var atoms_non_interface = master_atoms.Where(a => a.strand_interface_type == '_').ToList();
                        var atoms_non_interface_mpsa = atoms_non_interface.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_non_interface_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_non_interface.Select((a, i) => ((m == '*' || atoms_non_interface[i].dssp3_monomer == m || atoms_non_interface[i].dssp_monomer == m) && (atoms_non_interface[i].dssp3_monomer == atoms_non_interface_mpsa[i].predicted_ss_code || atoms_non_interface[i].dssp_monomer == atoms_non_interface_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_non_interface.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_non_interface_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_non_interface), a.ss, a.value, a.truth_total)));
                        var atoms_non_interface_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_non_interface.Select((a, i) => ((m == '*' || atoms_non_interface[i].stride3_monomer == m || atoms_non_interface[i].stride_monomer == m) && (atoms_non_interface[i].stride3_monomer == atoms_non_interface_mpsa[i].predicted_ss_code || atoms_non_interface[i].stride_monomer == atoms_non_interface_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_non_interface.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_non_interface_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_non_interface), a.ss, a.value, a.truth_total)));
                        var atoms_non_interface_mpsa_average = atoms_non_interface_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_non_interface_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_non_interface), a.ss, a.av, a.truth_total)));
                    }

                    {
                        // dimrophic
                        var atoms_dimorphic_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'S').ToList();
                        var atoms_dimorphic_strand_mpsa = atoms_dimorphic_strand.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_dimorphic_strand_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_dimorphic_strand.Select((a, i) => ((m == '*' || atoms_dimorphic_strand[i].dssp3_monomer == m || atoms_dimorphic_strand[i].dssp_monomer == m) && (atoms_dimorphic_strand[i].dssp3_monomer == atoms_dimorphic_strand_mpsa[i].predicted_ss_code || atoms_dimorphic_strand[i].dssp_monomer == atoms_dimorphic_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_dimorphic_strand.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_dimorphic_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_dimorphic_strand), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_dimorphic_strand.Select((a, i) => ((m == '*' || atoms_dimorphic_strand[i].stride3_monomer == m || atoms_dimorphic_strand[i].stride_monomer == m) && (atoms_dimorphic_strand[i].stride3_monomer == atoms_dimorphic_strand_mpsa[i].predicted_ss_code || atoms_dimorphic_strand[i].stride_monomer == atoms_dimorphic_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_dimorphic_strand.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_dimorphic_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_dimorphic_strand), a.ss, a.value, a.truth_total)));

                        var atoms_dimorphic_strand_mpsa_average = atoms_dimorphic_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_dimorphic_strand), a.ss, a.av, a.truth_total)));
                        var atoms_dimorphic_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'S') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'S') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_dimorphic_strand_flanking_mpsa = atoms_dimorphic_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_dimorphic_strand_flanking_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_dimorphic_strand_flanking.Select((a, i) => ((m == '*' || atoms_dimorphic_strand_flanking[i].dssp3_monomer == m || atoms_dimorphic_strand_flanking[i].dssp_monomer == m) && (atoms_dimorphic_strand_flanking[i].dssp3_monomer == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code || atoms_dimorphic_strand_flanking[i].dssp_monomer == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_dimorphic_strand_flanking.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_dimorphic_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_dimorphic_strand_flanking.Select((a, i) => ((m == '*' || atoms_dimorphic_strand_flanking[i].stride3_monomer == m || atoms_dimorphic_strand_flanking[i].stride_monomer == m) && (atoms_dimorphic_strand_flanking[i].stride3_monomer == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code || atoms_dimorphic_strand_flanking[i].stride_monomer == atoms_dimorphic_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_dimorphic_strand_flanking.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_dimorphic_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_mpsa_average = atoms_dimorphic_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_dimorphic_strand_flanking), a.ss, a.av, a.truth_total)));
                    }

                    {
                        // standard
                        var atoms_standard_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'M').ToList();
                        var atoms_standard_strand_mpsa = atoms_standard_strand.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_standard_strand_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_standard_strand.Select((a, i) => ((m == '*' || atoms_standard_strand[i].dssp3_monomer == m || atoms_standard_strand[i].dssp_monomer == m) && (atoms_standard_strand[i].dssp3_monomer == atoms_standard_strand_mpsa[i].predicted_ss_code || atoms_standard_strand[i].dssp_monomer == atoms_standard_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_standard_strand.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_standard_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_standard_strand), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_standard_strand.Select((a, i) => ((m == '*' || atoms_standard_strand[i].stride3_monomer == m || atoms_standard_strand[i].stride_monomer == m) && (atoms_standard_strand[i].stride3_monomer == atoms_standard_strand_mpsa[i].predicted_ss_code || atoms_standard_strand[i].stride_monomer == atoms_standard_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_standard_strand.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_standard_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_standard_strand), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_mpsa_average = atoms_standard_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_standard_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_standard_strand), a.ss, a.av, a.truth_total)));
                        var atoms_standard_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'M') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'M') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_standard_strand_flanking_mpsa = atoms_standard_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_standard_strand_flanking_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_standard_strand_flanking.Select((a, i) => ((m == '*' || atoms_standard_strand_flanking[i].dssp3_monomer == m || atoms_standard_strand_flanking[i].dssp_monomer == m) && (atoms_standard_strand_flanking[i].dssp3_monomer == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code || atoms_standard_strand_flanking[i].dssp_monomer == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_standard_strand_flanking.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_standard_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_standard_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_standard_strand_flanking_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_standard_strand_flanking.Select((a, i) => ((m == '*' || atoms_standard_strand_flanking[i].stride3_monomer == m || atoms_standard_strand_flanking[i].stride_monomer == m) && (atoms_standard_strand_flanking[i].stride3_monomer == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code || atoms_standard_strand_flanking[i].stride_monomer == atoms_standard_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_standard_strand_flanking.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_standard_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_standard_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_dimorphic_strand_flanking_mpsa_average = atoms_standard_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_dimorphic_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_standard_strand_flanking), a.ss, a.av, a.truth_total)));

                        //if (atoms_standard_strand.Count > 0)
                        //{
                        //    io_proxy.WriteLine(/*program.string_debug*/($@""));
                        //}
                    }


                    {
                        // hybrid
                        var atoms_hybrid_strand = master_atoms.Where((a, i) => a.strand_interface_type == 'H').ToList();
                        var atoms_hybrid_strand_mpsa = atoms_hybrid_strand.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_hybrid_strand_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_hybrid_strand.Select((a, i) => ((m == '*' || atoms_hybrid_strand[i].dssp3_monomer == m || atoms_hybrid_strand[i].dssp_monomer == m) && (atoms_hybrid_strand[i].dssp3_monomer == atoms_hybrid_strand_mpsa[i].predicted_ss_code || atoms_hybrid_strand[i].dssp_monomer == atoms_hybrid_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_hybrid_strand.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_hybrid_strand_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_hybrid_strand), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_hybrid_strand.Select((a, i) => ((m == '*' || atoms_hybrid_strand[i].stride3_monomer == m || atoms_hybrid_strand[i].stride_monomer == m) && (atoms_hybrid_strand[i].stride3_monomer == atoms_hybrid_strand_mpsa[i].predicted_ss_code || atoms_hybrid_strand[i].stride_monomer == atoms_hybrid_strand_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_hybrid_strand.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_hybrid_strand_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_hybrid_strand), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_mpsa_average = atoms_hybrid_strand_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_hybrid_strand_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_hybrid_strand), a.ss, a.av, a.truth_total)));
                        var atoms_hybrid_strand_flanking = master_atoms.Where((a, i) =>
                        {
                            // AA must not be strand interface (dimorphic/standard/hybrid)
                            if (a.strand_interface_type != '_') return false;

                            for (var j = 0; j <= 5; j++)
                            {
                                if (master_atoms.ElementAtOrDefault(i - j) != null && master_atoms[i - j].strand_interface_type == 'H') return true;
                                if (master_atoms.ElementAtOrDefault(i + j) != null && master_atoms[i + j].strand_interface_type == 'H') return true;
                            }

                            return false;
                        }).ToList();
                        var atoms_hybrid_strand_flanking_mpsa = atoms_hybrid_strand_flanking.SelectMany(a => a.mpsa_entries).Where(a => string.Equals(a.format, format, StringComparison.Ordinal)).Select(a => a.mpsa_entry).ToArray();
                        var atoms_hybrid_strand_flanking_q3_dssp = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_hybrid_strand_flanking.Select((a, i) => ((m == '*' || atoms_hybrid_strand_flanking[i].dssp3_monomer == m || atoms_hybrid_strand_flanking[i].dssp_monomer == m) && (atoms_hybrid_strand_flanking[i].dssp3_monomer == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code || atoms_hybrid_strand_flanking[i].dssp_monomer == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_hybrid_strand_flanking.Count(b => m == '*' || b.dssp3_monomer == m || b.dssp_monomer == m))).ToList();
                        atoms_hybrid_strand_flanking_q3_dssp.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"dssp"), format, nameof(atoms_hybrid_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_flanking_q3_stride = /*program.string_debug*/($@"HEC*").Select(m => (ss: m, value: (double)atoms_hybrid_strand_flanking.Select((a, i) => ((m == '*' || atoms_hybrid_strand_flanking[i].stride3_monomer == m || atoms_hybrid_strand_flanking[i].stride_monomer == m) && (atoms_hybrid_strand_flanking[i].stride3_monomer == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code || atoms_hybrid_strand_flanking[i].stride_monomer == atoms_hybrid_strand_flanking_mpsa[i].predicted_ss_code)) ? 1 : 0).Sum(), truth_total: (double)atoms_hybrid_strand_flanking.Count(b => m == '*' || b.stride3_monomer == m || b.stride_monomer == m))).ToList();
                        atoms_hybrid_strand_flanking_q3_stride.ForEach(a => q3.Add((pdb_id, chain_id, /*program.string_debug*/($@"stride"), format, nameof(atoms_hybrid_strand_flanking), a.ss, a.value, a.truth_total)));
                        var atoms_hybrid_strand_flanking_mpsa_average = atoms_hybrid_strand_flanking_mpsa.SelectMany(a => a.line_prob_values).GroupBy(a => a.ss).Select(a => (ss: a.Key, av: a.Select(b => b.value).Average(), truth_total: a.Count())).ToList();
                        atoms_hybrid_strand_flanking_mpsa_average.ForEach(a => av.Add((pdb_id, chain_id, format, nameof(atoms_hybrid_strand_flanking), a.ss, a.av, a.truth_total)));


                    }
                }


                //foreach (var s in mpsa_seqs.GroupBy(a => a.format))
                //{
                //    result1.Add(blank);

                //    format = s.Key;
                //    var x = s.ToList();

                //    // get the predicted SS sequence
                //    var predictor_seq = string.Join(/*program.string_debug*/($@""), x.Select(a => a.mpsa_entry.predicted_ss_code).ToArray());
                //    var predictor_seq_ss_types = predictor_seq.Union(x.First().mpsa_entry.ss_column_headers).Distinct().OrderBy(a => a).ToArray();
                //    result1.Add((filename, pdb, format, /*program.string_debug*/($@"seq_{format}", /*program.string_debug*/($@"seq_predicted", predictor_seq.Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));

                //    //predictor_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, /*program.string_debug*/($@"ss_binary", aa.ToString(CultureInfo.InvariantCulture), predictor_seq.Select(a => aa == a ? 1.0 : 0.0).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));

                //    predictor_seq_ss_types.ToList().ForEach(aa => result2.Add((filename, pdb, format, /*program.string_debug*/($@"ss_probability", aa.ToString(CultureInfo.InvariantCulture), x.Select(a => a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa).value).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray())));

                //    // distance between SS1 and SS1 for same amino acid on the same predictor - abs

                //    var distances = new List<(int index, double value)>();


                //    predictor_seq_ss_types.ToList().ForEach(aa1 => predictor_seq_ss_types.ToList().ForEach(aa2 =>
                //    {
                //        if (aa1 <= aa2) return;


                //        var d = x.Select((a, i) =>
                //        {
                //            var e = Math.Abs(a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa1).value - a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa2).value);

                //            distances.Add((i, e));

                //            return e;
                //        }).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray();


                //        result2.Add((filename, pdb, format, /*program.string_debug*/($@"ss_probability_distance_abs", /*program.string_debug*/($@"{aa1}_{aa2}", d));
                //    }));

                //    // 1. dimorphics interface average e%, h%, c%
                //    // 2. strand interface average e%, h%, c%
                //    // 3. protein average
                //    // 4. protein excl. strand interfaces average


                //    result2.Add((filename, pdb, format, /*program.string_debug*/($@"ss_probability_distance_abs", /*program.string_debug*/($@"average_all", distances.GroupBy(a => a.index).Select(a => a.Select(b => b.value).Average()).Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));


                // non-abs
                //predictor_seq_ss_types.ToList().ForEach(aa1 => predictor_seq_ss_types.ToList().ForEach(aa2 =>
                //{
                //    //if (aa1 <= aa2) return;
                //
                //    result2.Add((filename, pdb, format, /*program.string_debug*/($@"ss_probability_distance", /*program.string_debug*/($@"{aa1}_{aa2}",
                //        x.Select(a =>
                //                (a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa1).value -
                //                 a.mpsa_entry.line_prob_values.FirstOrDefault(b => b.ss == aa2).value))
                //            .Select(b => b.ToString(CultureInfo.InvariantCulture)).ToArray()));
                //}));


                //var score_dssp = 0d;
                //var score_stride = 0d;
                //for (var i = 0; i < dssp_seq.Length; i++)
                //{
                //    if (dssp3_seq[i] == predictor_seq[i] || dssp_seq[i] == predictor_seq[i]) score_dssp++;
                //    if (stride3_seq[i] == predictor_seq[i] || stride_seq[i] == predictor_seq[i]) score_stride++;

                //}

                //var dssp_pct = (double)score_dssp / (double)dssp_seq.Length;
                //var stride_pct = (double)score_stride / (double)dssp_seq.Length;

                //result1.Add((filename, pdb, format, /*program.string_debug*/($@"q3_dssp", /*program.string_debug*/($@"q3_dssp", new string[] { dssp_pct.ToString(CultureInfo.InvariantCulture) }));
                //result1.Add((filename, pdb, format, /*program.string_debug*/($@"q3_stride", /*program.string_debug*/($@"q3_stride", new string[] { stride_pct.ToString(CultureInfo.InvariantCulture) }));

                //    result2.Add(blank);
                //}

                //result1.Add(blank);
                //result2.Add(blank);

                //result1.GroupBy(a => a.filename).ToList().ForEach(a => program.AppendAllLines(Path.Combine(program.data_root_folder,"dssp_vs_mpsa\ss_predictor_comparison_{a.Key}.csv", a.Select(b => /*program.string_debug*/($@"{b.pdb},{b.format},{b.data_category1},{b.data_category2},{string.Join(/*program.string_debug*/($@","), b.data)}").ToList()));
                //result2.GroupBy(a => a.filename).ToList().ForEach(a => program.AppendAllLines(Path.Combine(program.data_root_folder,"dssp_vs_mpsa\ss_predictor_comparison_{a.Key}.csv", a.Select(b => /*program.string_debug*/($@"{b.pdb},{b.format},{b.data_category1},{b.data_category2},{string.Join(/*program.string_debug*/($@","), b.data)}").ToList()));

                var f_q3 = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"dssp_vs_mpsa"), /*program.string_debug*/($@"ss_predictor_comparison_{nameof(q3)}.csv"));
                var r_q3 = q3.Select(a => /*program.string_debug*/($@"{a.pdb_id},{a.chain_id},{a.truth_format},{a.predictor_format},{a.aa_subset},{a.ss},{a.q3_value},{a.truth_total}")).ToList();

                io_proxy.AppendAllLines(f_q3, r_q3, nameof(atom), nameof(compare_dssp_to_mpsa));

                var f_av = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"dssp_vs_mpsa"), /*program.string_debug*/($@"ss_predictor_comparison_{nameof(av)}.csv"));
                var r_av = av.Select(a => /*program.string_debug*/($@"{a.pdb_id},{a.chain_id},{a.predictor_format},{a.aa_subset},{a.ss},{a.average},{a.truth_total}")).ToList();

                io_proxy.AppendAllLines(f_av, r_av, nameof(atom), nameof(compare_dssp_to_mpsa));
            }
        }



        internal static List<(string format, string prediction)> get_dssp_and_mpsa_subsequences(List<atom> atoms, enum_get_dssp_and_mpsa_subsequences_params get_dssp_and_mpsa_subsequences_params = enum_get_dssp_and_mpsa_subsequences_params.none)
        {
            var result = new List<(string format, string prediction)>();

            //result.Add(/*program.string_debug*/($@"DSSP vs MPSA: " + pdb_id + /*program.string_debug*/($@" ") + chain_id);

            //atoms = atoms.Where(a => a.pdb_id == pdb_id && a.chain_id == chain_id).ToList();

            var master_atoms = atom.select_amino_acid_master_atoms(null, atoms);

            if (master_atoms == null || master_atoms.Count == 0) return result;

            var aa_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.amino_acid).ToList());

            var monomer_dssp_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp_monomer).ToList());
            var monomer_stride_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride_monomer).ToList());
            var monomer_dssp3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp3_monomer).ToList());
            var monomer_stride3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride3_monomer).ToList());

            var multimer_dssp_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp_multimer).ToList());
            var multimer_stride_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride_multimer).ToList());
            var multimer_dssp3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.dssp3_multimer).ToList());
            var multimer_stride3_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a => a.stride3_multimer).ToList());

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.aa_seq)) { result.Add((nameof(aa_seq), aa_seq)); }

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.monomer_dssp_seq)) { result.Add((nameof(monomer_dssp_seq), monomer_dssp_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.monomer_stride_seq)) { result.Add((nameof(monomer_stride_seq), monomer_stride_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.monomer_dssp3_seq)) { result.Add((nameof(monomer_dssp3_seq), monomer_dssp3_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.monomer_stride3_seq)) { result.Add((nameof(monomer_stride3_seq), monomer_stride3_seq)); }

            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.multimer_dssp_seq)) { result.Add((nameof(multimer_dssp_seq), multimer_dssp_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.multimer_stride_seq)) { result.Add((nameof(multimer_stride_seq), multimer_stride_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.multimer_dssp3_seq)) { result.Add((nameof(multimer_dssp3_seq), multimer_dssp3_seq)); }
            if (get_dssp_and_mpsa_subsequences_params.HasFlag(enum_get_dssp_and_mpsa_subsequences_params.multimer_stride3_seq)) { result.Add((nameof(multimer_stride3_seq), multimer_stride3_seq)); }


            var mpsa_seqs = master_atoms.SelectMany(a => a?.mpsa_entries ?? new List<(string format, info_mpsa_reader_line_entry mpsa_entry)>()).ToList();

            foreach (var s in mpsa_seqs.GroupBy(a => a.format))
            {
                var format = s.Key;
                var x = s.ToList();
                var mpsa_seq = string.Join(/*program.string_debug*/($@""), x.Select(a => a.mpsa_entry.predicted_ss_code).ToList());

                result.Add((format, mpsa_seq));
            }

            var con_ss_seq = string.Join(/*program.string_debug*/($@""), master_atoms.Select(a =>
            {
                var y = a?.mpsa_entries?.SelectMany(b => b.mpsa_entry?.line_prob_values ?? new List<(char ss, char amino_acid, double value)>()).ToList() ?? new List<(char ss, char amino_acid, double value)>();

                var prob_h = y.Where(b => b.ss == 'H').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_e = y.Where(b => b.ss == 'E').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_c = y.Where(b => b.ss == 'C').Select(b => b.value).DefaultIfEmpty(0).Average();
                var prob_t = y.Where(b => b.ss == 'T').Select(b => b.value).DefaultIfEmpty(0).Average();
                //var prob_f = y.Where(b => b.ss == 'F').Select(b => b.value).DefaultIfEmpty(0).Average();

                var values = new List<(char ss, double value)> { ('C', prob_c), ('E', prob_e), ('H', prob_h), ('T', prob_t) };
                //values.Add(('F', prob_f));
                values = values.OrderByDescending(c => c.value).ToList();

                var all_zero = values.All(d => d.value == 0);

                return all_zero ? 'C' : values.First().ss;
            }).ToList());

            result.Add((/*program.string_debug*/($@"consensus"), con_ss_seq));

            return result;
        }


        internal static List<atom> select_amino_acid_master_atoms(string pdb_id, List<atom> atoms, string atom_type = null)
        {
            // if null, default to CA or first available in N, C, O, CB.

            if (atoms == null || atoms.Count == 0)
            {
                return atoms;
            }

            if (!string.IsNullOrWhiteSpace(pdb_id))
            {
                atoms = atoms.Where(a => /*pdb_id == null || */string.Equals(pdb_id, a.pdb_id, StringComparison.OrdinalIgnoreCase)).Distinct().ToList();
            }

            var master_atoms = atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code)
            //    .ToList();
            //atoms = atoms

            .GroupBy(atom => (atom.pdb_id, atom.chain_id, atom.residue_index /*, atom.i_code*/))
                .Select(grouped_atoms =>
                {
                    if (!string.IsNullOrWhiteSpace(atom_type))
                    {
                        var rq = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, atom_type, StringComparison.Ordinal));

                        return rq; // return rq, even if null, as if does not exist, do not want replacement/substitute atom
                    }

                    if (grouped_atoms.Count() == 1)
                    {
                        return grouped_atoms.First();
                    }

                    var ca = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"CA"), StringComparison.Ordinal));
                    if (ca != null) return ca;

                    var n = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"N"), StringComparison.Ordinal));
                    if (n != null) return n;

                    var c = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"C"), StringComparison.Ordinal));
                    if (c != null) return c;

                    var o = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"O"), StringComparison.Ordinal));
                    if (o != null) return o;

                    var cb = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"CB"), StringComparison.Ordinal));
                    if (cb != null) return cb;

                    return grouped_atoms.FirstOrDefault(atom => atom != null);
                }).Where(a => a != null).ToList();

            //master_atoms = master_atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

            return master_atoms;
        }

        internal static List<atom> select_amino_acid_master_atoms_ignore_chain(string pdb_id, List<atom> atoms)
        {
            atoms = atoms.Where(a => pdb_id == null || string.Equals(pdb_id, a.pdb_id, StringComparison.OrdinalIgnoreCase)).Distinct().ToList();

            atoms = atoms.OrderBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();

            atoms = atoms.GroupBy(atom => (atom.pdb_id, atom.residue_index /*, atom.i_code*/)).Select(grouped_atoms =>
            {
                var ca = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"CA"), StringComparison.Ordinal));
                if (ca != null) return ca;

                var n = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"N"), StringComparison.Ordinal));
                if (n != null) return n;

                var c = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"C"), StringComparison.Ordinal));
                if (c != null) return c;

                var o = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"O"), StringComparison.Ordinal));
                if (o != null) return o;

                var cb = grouped_atoms.FirstOrDefault(b => string.Equals(b.atom_type, /*program.string_debug*/($@"CB"), StringComparison.Ordinal));
                if (cb != null) return cb;

                return grouped_atoms.FirstOrDefault(atom => atom != null);
            }).ToList();

            return atoms;
        }

        internal static List<(string pdb_id, char chain_id, List<atom> atom_list)> group_atoms_by_pdb_id_and_chain_id(string pdb_id, List<atom> atoms)
        {
            atoms = atoms.Distinct().ToList();
            atoms = atoms.OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code) /*.ThenBy(a => a.SerialIndex)*/.ToList();

            var seq = atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.OrdinalIgnoreCase)).GroupBy(a => (pdb_id: a.pdb_id.ToUpperInvariant(), chain_id: a.chain_id)).Select(a => (pdb_id: a.Key.pdb_id, chain_id: a.Key.chain_id, atom_list: a.ToList())).ToList();
            return seq;
        }

        internal static List<(string pdb_id, char chain_id, List<atom> atom_list)> group_master_atoms_by_pdb_id_and_chain_id(string pdb_id, List<atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);


            atoms = select_amino_acid_master_atoms(pdb_id, atoms);
            var seq = atoms.GroupBy(a => (pdb_id: a.pdb_id.ToUpperInvariant(), chain_id: a.chain_id)).Select(a => (pdb_id: a.Key.pdb_id, chain_id: a.Key.chain_id, atom_list: a.ToList())).ToList();
            return seq;
        }

        internal static List<(string pdb_id, char chain_id, string aa_sequence)> amino_acid_sequence(string pdb_id, List<atom> atoms)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var seq = group_master_atoms_by_pdb_id_and_chain_id(pdb_id, atoms).Select(a => (pdb_id: a.pdb_id, chain_id: a.chain_id, aa_sequence: string.Join(/*program.string_debug*/($@""), a.atom_list.Select(b => b.amino_acid).ToList()))).ToList();
            return seq;
        }

        internal static List<(string pdb_id, char chain_id, string ss_sequence)> ss_sequence(string pdb_id, List<atom> atoms, enum_structure_oligomisation structure_oligomisation, enum_ss_type ss_type = enum_ss_type.DSSP)
        {
            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);


            return group_master_atoms_by_pdb_id_and_chain_id(pdb_id, atoms).Select(a => (pdb_id: a.pdb_id, chain_id: a.chain_id, ss_sequence: string.Join(/*program.string_debug*/($@""), a.atom_list.Select(b =>
            {
                if (ss_type == enum_ss_type.DSSP) return structure_oligomisation == enum_structure_oligomisation.multimer ? b.dssp_multimer : b.dssp_monomer;
                if (ss_type == enum_ss_type.DSSP3) return structure_oligomisation == enum_structure_oligomisation.multimer ? b.dssp3_multimer : b.dssp3_monomer;

                if (ss_type == enum_ss_type.STRIDE) return structure_oligomisation == enum_structure_oligomisation.multimer ? b.stride_multimer : b.stride_monomer;
                if (ss_type == enum_ss_type.STRIDE3) return structure_oligomisation == enum_structure_oligomisation.multimer ? b.stride3_multimer : b.stride3_monomer;
                return ' ';
            }).ToList()))).ToList();
        }

        private static void load_mpsa_sec_struct_predictions(string pdb_id, List<atom> pdb_model_atoms)
        {
            //this.mpsa_entries = new mpsa_reader(/*program.string_debug*/($@""));

            var pdb_id_simple = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            var master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);


            //var lookup_sequence_table = io_proxy.ReadAllLines(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"betastrands_dataset_sequences.txt"), nameof(atom), nameof(load_mpsa_sec_struct_predictions)).ToList();


            var sequences = atom.amino_acid_sequence(pdb_id, master_atoms);

            foreach (var seq in sequences)
            {
                var chain_id = seq.chain_id;

                var chain_master_atoms = master_atoms.Where(a => a.chain_id == chain_id).ToList();
                //var sequence = seq.aa_sequence;

                //var line_index_match = lookup_sequence_table.IndexOf(sequence);

                //if (line_index_match < 0) throw new Exception();

                //if (pdb_id_simple == /*program.string_debug*/($@"1EGP")
                //{
                //}

                var format_list = info_mpsa_reader.secondary_structure_codes.Select(a => a.format).ToList();

                //var mpsa_files_wc = Path.Combine(program.data_root_folder,"ss_meta\{line_index_match}.*";
                //var mpsa_files = Directory.GetFiles(Path.GetDirectoryName(mpsa_files_wc), Path.GetFileName(mpsa_files_wc));

                var mpsa_files = format_list.Select(a => Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ss_meta"), /*program.string_debug*/($@"{pdb_id_simple}{chain_id}.{a}"))).ToList();

                foreach (var mpsa_file in mpsa_files)
                {
                    var reader = new info_mpsa_reader(mpsa_file, chain_master_atoms.Select(a => a.amino_acid).ToList());


                    for (var index = 0; index < chain_master_atoms.Count; index++)
                    {
                        var master_atom = chain_master_atoms[index];

                        if (master_atom.mpsa_entries == null) master_atom.mpsa_entries = new List<(string format, info_mpsa_reader_line_entry)>();


                        var mpsa_matrix_line_entry = reader?.mpsa_matrix?.First /*OrDefault*/(a => a.index == index);

                        //if (mpsa_matrix == null || mpsa_matrix.line_prob_values==null || mpsa_matrix.line_prob_values.Count==0)
                        //{
                        //    mpsa_matrix = new mpsa_reader.mpsa_line_entry()
                        //    {   reader = reader,
                        //        index = index,
                        //        predicted_ss_code = ' ',
                        //        amino_acid = master_atom.amino_acid,
                        //        line_prob_values = reader.ss_column_headers.Select(a=>(a, master_atom.amino_acid, 0d)).ToList(),
                        //        ss_column_headers = reader.ss_column_headers
                        //    };
                        //}

                        master_atom.mpsa_entries.Add((reader.format, mpsa_matrix_line_entry));

                        foreach (var atom in master_atom.amino_acid_atoms)
                        {
                            atom.mpsa_entries = master_atom.mpsa_entries;
                        }
                    }

                    //}
                }

                // build consensus
            }
        }


        //internal static List<string> pssm_database_names = new List<string>();

        internal static void load_dna_binding_stackdppred(string pdb_id, List<atom> pdb_model_atoms)
        {
            if (pdb_model_atoms == null || pdb_model_atoms.Count == 0)
            {
                return;
            }

            var files = Directory.GetFiles(Path.Combine(program.data_root_folder, /*program.string_debug*/($@"stackdppred_results")), /*program.string_debug*/($@"{pdb_id}.txt"), SearchOption.AllDirectories);

            foreach (var file in files)
            {
                var data = io_proxy.ReadAllLines(file, nameof(atom), nameof(load_dna_binding_stackdppred)).First();
                var data2 = string.Join(/*program.string_debug*/($@""), data.Where(a => !"{}'':,".Contains(a, StringComparison.Ordinal)).ToList()).Split().ToList();

                var non_binding_prob = double.Parse(data2[data2.IndexOf(/*program.string_debug*/($@"non-binding_prob")) + 1], NumberStyles.Float, NumberFormatInfo.InvariantInfo);
                var binding_prob = double.Parse(data2[data2.IndexOf(/*program.string_debug*/($@"binding_prob")) + 1], NumberStyles.Float, NumberFormatInfo.InvariantInfo);

                var db = Path.GetDirectoryName(file)?.Split(new char[] { '\\', '/' }, StringSplitOptions.RemoveEmptyEntries).Last().Split('_').Last();

                if (string.IsNullOrWhiteSpace(db)) continue;

                foreach (var atom in pdb_model_atoms)
                {
                    if (string.Equals(db, /*program.string_debug*/($@"nr"), StringComparison.OrdinalIgnoreCase)) atom.chain_dna_binding_prob_nr = binding_prob;
                    else if (string.Equals(db, /*program.string_debug*/($@"swissprot"), StringComparison.OrdinalIgnoreCase)) atom.chain_dna_binding_prob_swissprot = binding_prob;
                    else if (string.Equals(db, /*program.string_debug*/($@"uniref90"), StringComparison.OrdinalIgnoreCase)) atom.chain_dna_binding_prob_uniref90 = binding_prob;
                }

                //io_proxy.WriteLine();

            }
        }

        internal static void load_pssm(string pdb_id, List<atom> pdb_model_atoms)
        {
            var pdb_id_simple = Path.GetFileNameWithoutExtension(pdb_id).Substring(0, 4);

            var master_atoms = atom.select_amino_acid_master_atoms(pdb_id, pdb_model_atoms);

            //var lookup_sequence_table = program.ReadAllLines(Path.Combine(program.data_root_folder,"betastrands_dataset_sequences.txt").ToList();

            var sequences = atom.amino_acid_sequence(pdb_id, master_atoms);

            foreach (var seq in sequences)
            {
                var chain_id = seq.chain_id;

                var chain_master_atoms = master_atoms.Where(a => a.chain_id == chain_id).ToList();
                var sequence = seq.aa_sequence;

                //var line_index_match = lookup_sequence_table.IndexOf(sequence);

                //if (line_index_match < 0) throw new Exception();

                var pssm_folders = Directory.GetDirectories(Path.Combine(program.data_root_folder), /*program.string_debug*/($@"blast_pssm_*"));
                var pssm_database_names = pssm_folders.Select(a => a.Split(new char[] { '\\', '/' }, StringSplitOptions.RemoveEmptyEntries).Last()).ToList();


                //var pssm_files = new string[]
                //    {
                //        Path.Combine(program.data_root_folder,"blast_pssm_swissprot_local\{line_index_match}.pssm",
                //        Path.Combine(program.data_root_folder,"blast_pssm_nr_remote\{line_index_match}.pssm"
                //    };

                var pssm_files = pssm_folders.Select(a => Path.Combine(a, /*program.string_debug*/($@"{pdb_id_simple}{chain_id}.pssm"))).ToArray();

                foreach (var pssm_file in pssm_files)
                {
                    var pssm_database_name = Path.GetDirectoryName(pssm_file).Split(new char[] { '\\', '/' }, StringSplitOptions.RemoveEmptyEntries).Last();

                    //if (File.Exists(pssm_file) && new FileInfo(pssm_file).Length > 0)
                    {
                        var pssm_matrix_unnormalised = info_blast_pssm.load_psi_blast_pssm(pssm_file);
                        var pssm_matrix_normalised = info_blast_pssm.normalise_pssm(pssm_matrix_unnormalised);

                        for (var index = 0; index < chain_master_atoms.Count; index++)
                        {
                            //var amino_acid_pssm = pssm_match.First(a => a.matrix_row_index == index);
                            var amino_acid_pssm_unnormalised = pssm_matrix_unnormalised.Where(a => a.matrix_row_index == index).ToList();
                            var amino_acid_pssm_normalised = pssm_matrix_normalised.Where(a => a.matrix_row_index == index).ToList();

                            if (amino_acid_pssm_unnormalised.Count != 20)
                            {
                                //throw new Exception();
                                amino_acid_pssm_unnormalised = new List<info_blast_pssm_entry>();
                                amino_acid_pssm_normalised = new List<info_blast_pssm_entry>();
                            }

                            var master_atom = chain_master_atoms[index];

                            if (master_atom.amino_acid_pssm_unnormalised == null) master_atom.amino_acid_pssm_unnormalised = new List<(string database, List<info_blast_pssm_entry> pssm_entries)>();
                            if (master_atom.amino_acid_pssm_normalised == null) master_atom.amino_acid_pssm_normalised = new List<(string database, List<info_blast_pssm_entry> pssm_entries)>();

                            master_atom.amino_acid_pssm_unnormalised.Add((pssm_database_name, amino_acid_pssm_unnormalised));
                            master_atom.amino_acid_pssm_normalised.Add((pssm_database_name, amino_acid_pssm_normalised));

                            foreach (var a in master_atom.amino_acid_atoms)
                            {
                                a.amino_acid_pssm_unnormalised = master_atom.amino_acid_pssm_unnormalised;
                                a.amino_acid_pssm_normalised = master_atom.amino_acid_pssm_normalised;
                            }
                        }
                    }
                }
            }
        }

        internal static void load_ring(string pdb_id, List<atom> atoms)
        {
            //var edge_files = Directory.GetFiles(@"C:\phd\ring_pdb\pdb\", /*program.string_debug*/($@"*.edges");
            //var node_files = Directory.GetFiles(@"C:\phd\ring_pdb\pdb\", /*program.string_debug*/($@"*.nodes");

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var ring_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ring_pdb"), /*program.string_debug*/($@"pdb"));
            // Note: careful from the added '_Repair' in the filename

            var chains = atoms.Select(a => a.chain_id).Distinct().ToList();

            foreach (var chain in chains)
            {
                {
                    var ring_monomer_edge_filename = Path.Combine(ring_folder, /*program.string_debug*/($@"{pdb_id}{chain}_Repair.pdb.edges"));

                    if (!File.Exists(ring_monomer_edge_filename) || new FileInfo(ring_monomer_edge_filename).Length == 0)
                    {
                        var missing_edges_file = Path.Combine(ring_folder, /*program.string_debug*/($@"missing_edges.txt"));
                        var missing_edges_text = /*program.string_debug*/($@"Warning: Ring Edge File is missing or empty: {ring_monomer_edge_filename}");

                        io_proxy.AppendAllLines(missing_edges_file, new string[] { missing_edges_text }, nameof(atom), nameof(load_ring));
                        io_proxy.WriteLine(missing_edges_text, nameof(atom), nameof(load_ring));
                    }

                    var ring_monomer_edges = info_ring_edge.load(ring_monomer_edge_filename);

                    foreach (var ring_monomer_edge in ring_monomer_edges)
                    {
                        var edge_atoms = atoms.Where(a => (
                                ring_monomer_edge.NodeId1.chain == a.chain_id &&
                                ring_monomer_edge.NodeId1.res_id == a.residue_index &&
                                ring_monomer_edge.NodeId1.icode == a.i_code &&
                                ring_monomer_edge.NodeId1.amino_acid1 == a.amino_acid
                                )
                        //||
                        //(ring_monomer_edge.NodeId2.chain == a.chain_id &&
                        // ring_monomer_edge.NodeId2.res_id == a.residue_index &&
                        // ring_monomer_edge.NodeId2.icode == a.i_code &&
                        // ring_monomer_edge.NodeId2.amino_acid1 == a.amino_acid)
                        ).ToList();

                        //if (edge_atoms == null || edge_atoms.Count == 0)
                        //{
                        //    io_proxy.WriteLine();
                        //}

                        foreach (var a in edge_atoms)
                        {
                            if (a.monomer_ring_edges == null) a.monomer_ring_edges = new List<info_ring_edge>();
                            a.monomer_ring_edges.Add(ring_monomer_edge);
                        }
                    }
                }

                {
                    var ring_monomer_node_filename = Path.Combine(ring_folder, /*program.string_debug*/($@"{pdb_id}{chain}_Repair.pdb.nodes"));

                    if (!File.Exists(ring_monomer_node_filename) || new FileInfo(ring_monomer_node_filename).Length == 0)
                    {
                        //io_proxy.WriteLine(/*program.string_debug*/($@"Warning: Ring Node File is missing or empty: {ring_monomer_node_filename}");

                        var missing_nodes_file = Path.Combine(ring_folder, /*program.string_debug*/($@"missing_nodes.txt"));
                        var missing_nodes_text = /*program.string_debug*/($@"Warning: Ring Node File is missing or empty: {ring_monomer_node_filename}");

                        io_proxy.AppendAllLines(missing_nodes_file, new string[] { missing_nodes_text }, nameof(atom), nameof(load_ring));
                        io_proxy.WriteLine(missing_nodes_text, nameof(atom), nameof(load_ring));

                    }

                    var ring_monomer_nodes = info_ring_node.load(ring_monomer_node_filename);

                    foreach (var ring_monomer_node in ring_monomer_nodes)
                    {
                        var node_atoms = atoms.Where(a => ring_monomer_node.Chain == a.chain_id && ring_monomer_node.Position == a.residue_index && ring_monomer_node.Residue1 == a.amino_acid).ToList();

                        //if (node_atoms == null || node_atoms.Count == 0)
                        //{
                        //    io_proxy.WriteLine();
                        //}

                        foreach (var a in node_atoms)
                        {
                            if (a.monomer_ring_nodes == null) a.monomer_ring_nodes = new List<info_ring_node>();
                            a.monomer_ring_nodes.Add(ring_monomer_node);
                        }
                    }
                }
            }
        }

        internal static void load_dssp(string pdb_id, List<atom> atoms)
        {
            if (atoms == null || atoms.Count == 0)
            {
                return;
            }

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var multimer_dssp_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ss_dssp"), /*program.string_debug*/($@"{pdb_id}.dssp"));

            var multimer_dssp_list = File.Exists(multimer_dssp_list_file) ? info_dssp.Load(multimer_dssp_list_file).ToList() : new List<info_dssp_item>();

            var chain_ids = atoms.Select(a => a.chain_id).Distinct().ToList();


            var monomer_dssp_list = chain_ids.SelectMany(chain_id =>
            {
                var monomer_dssp_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ss_dssp"), /*program.string_debug*/($@"{pdb_id}{chain_id}.dssp"));
                if (File.Exists(monomer_dssp_list_file))
                {
                    return info_dssp.Load(monomer_dssp_list_file).ToList();
                }
                else
                {
                    return new List<info_dssp_item>();
                }
            }).ToList();

            //const char default_dssp_code = ' ';

            foreach (var atom in atoms)
            {
                {
                    var multimer_dssp = multimer_dssp_list.FirstOrDefault(a => (a.Chain == atom.chain_id) && (a.PdbResidueSequenceIndex.Trim() == atom.residue_index.ToString(CultureInfo.InvariantCulture).Trim()) && (a.iCode == atom.i_code));

                    atom.dssp_multimer = multimer_dssp?.SecondaryStructure ?? info_dssp_item.dssp_default_secondary_structure;
                    atom.dssp3_multimer = secondary_structure_state_reduction(atom.dssp_multimer);
                }
                {
                    var monomer_dssp = monomer_dssp_list.FirstOrDefault(a => (a.Chain == atom.chain_id) && (a.PdbResidueSequenceIndex.Trim() == atom.residue_index.ToString(CultureInfo.InvariantCulture).Trim()) && (a.iCode == atom.i_code));

                    atom.dssp_monomer = monomer_dssp?.SecondaryStructure ?? info_dssp_item.dssp_default_secondary_structure;
                    atom.dssp3_monomer = secondary_structure_state_reduction(atom.dssp_monomer);
                }
            }
        }


        internal static void load_stride(string pdb_id, List<atom> atoms)
        {
            if (atoms == null || atoms.Count == 0)
            {
                return;
            }

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            var stride_multimer_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ss_stride"), /*program.string_debug*/($@"{pdb_id}.stride"));

            var stride_multimer_list = File.Exists(stride_multimer_list_file) ? info_stride.Load(stride_multimer_list_file).Where(a => a.GetType() == typeof(info_stride.Stride_DetailedSecondaryStructureAssignments)).Select(a => (info_stride.Stride_DetailedSecondaryStructureAssignments)a).ToList() : new List<info_stride.Stride_DetailedSecondaryStructureAssignments>();

            var chain_ids = atoms.Select(a => a.chain_id).Distinct().ToList();

            var stride_monomer_list = chain_ids.SelectMany(chain_id =>
            {
                var stride_monomer_list_file = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"ss_stride"), /*program.string_debug*/($@"{pdb_id}{chain_id}.stride"));

                if (File.Exists(stride_monomer_list_file))
                {
                    return info_stride.Load(stride_monomer_list_file);
                }
                else
                {
                    return new List<info_stride.stride_record>();
                }
            }).Where(a => a.GetType() == typeof(info_stride.Stride_DetailedSecondaryStructureAssignments)).Select(a => (info_stride.Stride_DetailedSecondaryStructureAssignments)a).ToList();

            foreach (var atom in atoms)
            {
                {
                    var stride_multimer = stride_multimer_list.FirstOrDefault(a => (string.Equals(a.ProteinChainIdentifier.Trim(), atom.chain_id.ToString(CultureInfo.InvariantCulture).Trim(), StringComparison.Ordinal)) && (string.Equals(a.PdbResidueNumber.Trim(), atom.residue_index.ToString(CultureInfo.InvariantCulture).Trim(), StringComparison.Ordinal)));

                    atom.stride_multimer = stride_multimer?.OneLetterSecondaryStructureCode?[0] ?? info_stride.Stride_DetailedSecondaryStructureAssignments.stride_default_secondary_structure;
                    atom.stride3_multimer = secondary_structure_state_reduction(atom.stride_multimer);
                }
                {
                    var stride_monomer = stride_monomer_list.FirstOrDefault(a => (string.Equals(a.ProteinChainIdentifier.Trim(), atom.chain_id.ToString(CultureInfo.InvariantCulture).Trim(), StringComparison.Ordinal)) && (string.Equals(a.PdbResidueNumber.Trim(), atom.residue_index.ToString(CultureInfo.InvariantCulture).Trim(), StringComparison.Ordinal)));

                    atom.stride_monomer = stride_monomer?.OneLetterSecondaryStructureCode?[0] ?? info_stride.Stride_DetailedSecondaryStructureAssignments.stride_default_secondary_structure;
                    atom.stride3_monomer = secondary_structure_state_reduction(atom.stride_monomer);
                }
            }
        }

        internal static void load_rsa(string pdb_id, List<atom> atoms)
        {
            if (atoms == null || atoms.Count == 0)
            {
                return;
            }

            pdb_id = Path.GetFileNameWithoutExtension(pdb_id);

            //atoms.Select(a => a.PdbId).ToList();

            var rsa_folder = Path.Combine(program.data_root_folder, /*program.string_debug*/($@"sasa"));

            if (!Directory.Exists(rsa_folder)) return;

            var rsa_file_wildcard = /*program.string_debug*/($@"{pdb_id}*.?-rsa-sasa");

            var rsa_files = Directory.GetFiles(rsa_folder, rsa_file_wildcard).ToList();
            var rsa_list = info_solvent_access.load(rsa_files);

            foreach (var atom in atoms)
            {
                atom.RSA_L = rsa_list.FirstOrDefault(a => (string.Equals(a.pdb_id, atom.pdb_id, StringComparison.OrdinalIgnoreCase)) && (a.chain_id == atom.chain_id) && (a.algo == 'L') && (a.res_num == atom.residue_index) && (a.amino_acid == atom.amino_acid));
                atom.RSA_S = rsa_list.FirstOrDefault(a => (string.Equals(a.pdb_id, atom.pdb_id, StringComparison.OrdinalIgnoreCase)) && (a.chain_id == atom.chain_id) && (a.algo == 'S') && (a.res_num == atom.residue_index) && (a.amino_acid == atom.amino_acid));
            }
        }

        internal static double Distance3D(atom atom1, atom atom2)
        {
            //if (atom1 == null)
            //{
            //    throw new ArgumentNullException(nameof(atom1));
            //}

            //if (atom2 == null)
            //{
            //    throw new ArgumentNullException(nameof(atom2));
            //}

            return Distance3D(atom1.X, atom1.Y, atom1.Z, atom2.X, atom2.Y, atom2.Z);
        }

        internal static double Distance2D(double x0, double y0, double x1, double y1)
        {
            return (double)Math.Sqrt((double)(((x1 - x0) * (x1 - x0)) + ((y1 - y0) * (y1 - y0))));
        }


        internal static double Distance3D(double x0, double y0, double z0, double x1, double y1, double z1)
        {
            return (double)Math.Sqrt((double)(((x1 - x0) * (x1 - x0)) + ((y1 - y0) * (y1 - y0)) + ((z1 - z0) * (z1 - z0))));
        }


        internal static double measure_atomic_linear_distance(List<atom> atoms)
        {
            //atoms = select_amino_acid_master_atoms(null, atoms);

            var atom1 = atoms.First();
            var atom2 = atoms.Last();

            var distance = Distance3D(atom1.X, atom1.Y, atom1.Z, atom2.X, atom2.Y, atom2.Z);

            return distance;
        }

        internal static double measure_atomic_curve_distance(List<atom> atoms)
        {
            // atoms must be in correct order
            // measures the distance in angstrom of the curve
            //atoms = select_amino_acid_master_atoms(null, atoms);

            var curve_distances = new double[atoms.Count - 1];
            for (var index = 0; index < atoms.Count - 1; index++)
            {
                var atom1 = atoms[index];
                var atom2 = atoms[index + 1];

                curve_distances[index] = Distance3D(atom1.X, atom1.Y, atom1.Z, atom2.X, atom2.Y, atom2.Z);
            }

            return curve_distances.Sum();
        }

        internal static (double displacement, double distance_of_curve, double tortuosity1) measure_tortuosity1(List<atom> atoms)
        {
            if (atoms == null || atoms.Count <= 1)
            {
                return (0, 0, 0);
            }


            //atoms = select_amino_acid_master_atoms(null, atoms); // note: atoms parameter should already be passed as master atoms

            var distance_of_curve = measure_atomic_curve_distance(atoms);

            //var distance_of_curve = distance_of_curve_stats.sum;

            var displacement = measure_atomic_linear_distance(atoms);

            var tortuosity1 = displacement == 0 ? 0d : (distance_of_curve / displacement);

            return (displacement, distance_of_curve, tortuosity1);
        }

        internal static List<T> list_subsection_from_obj_to_obj<T>(List<T> list, T first_list_object, T last_list_object)
        {
            if (list == null)
            {
                return null;
            }

            var index1 = list.IndexOf(first_list_object);
            var index2 = list.IndexOf(last_list_object);

            return list.Skip(index1 < index2 ? index1 : index2).Take(Math.Abs(index2 - index1) + 1).ToList();
        }

        //internal static (descriptive_stats displacement_stat_values, descriptive_stats curve_stat_values, descriptive_stats tortuosity_stat_values) measure_tortuosity2(List<atom> atoms)
        internal static (double[] displacement_stat_values, double[] curve_stat_values, double[] tortuosity_stat_values) measure_tortuosity2(List<atom> atoms)
        {
            //if (atoms == null || atoms.Count <= 1)
            //{
            //    return (null, null, null);
            //}

            //atoms = select_amino_acid_master_atoms(null, atoms);



            var distances = Array.Empty<(double curve_distance, double displacement_distance)>();

            if (atoms != null && atoms.Count >= 2)
            {
                var pairs = atoms.SelectMany((a, x) => atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).ToList()).ToList();
                distances = pairs.Select(a =>
                    (
                        curve_distance: measure_atomic_curve_distance(list_subsection_from_obj_to_obj(atoms, a.atom1, a.atom2)),
                        displacement_distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z)
                        )).ToArray();
            }

            // linear_distance / linear_list / linear_stat_values --> list of distances between all pairs of atoms, then averaged
            // curve_distance / 
            var displacement_distance_list = distances.Select(a => a.displacement_distance).OrderBy(a=>a).ToArray(); // list of LineOfSight distances between all points (CA atoms)
            //var curve_stat_list = distances.Select(a => a.curve_distance).ToArray(); // list of all curve distances between all pairs of points (CA atoms) 
            var curve_distance_list = distances.Select(a => a.curve_distance).OrderBy(a=>a).ToArray(); //curve_stat_list.Select(a => a.sum).ToArray();

            var tortuosity_list = distances.Select(a => a.displacement_distance == 0 ? 0d : a.curve_distance / a.displacement_distance).OrderBy(a=>a).ToArray();

            return (displacement_distance_list, curve_distance_list, tortuosity_list);

          

            //return (displacement_stat_values, curve_stat_values, tortuosity_stat_values);
        }
        /*
        internal static void find_intramolecular_contacts(string pdb_id, List<atom> atoms, double? max_dist = null)
        {
            if (atoms == null)
            {
                throw new ArgumentNullException(nameof(atoms));
            }

            const bool intermolecular = false;
            const bool intramolecular = true;
            find_contacts(pdb_id, atoms, intermolecular: intermolecular, intramolecular: intramolecular, max_dist: max_dist);
        }

        internal static void find_intermolecular_contacts(string pdb_id, List<atom> atoms, double? max_dist = null)
        {
            if (atoms == null)
            {
                throw new ArgumentNullException(nameof(atoms));
            }

            const bool intermolecular = true;
            const bool intramolecular = false;
            find_contacts(pdb_id, atoms, intermolecular: intermolecular, intramolecular: intramolecular, max_dist: max_dist);
        }

        
        private static void find_contacts(string pdb_id, List<atom> atoms, bool intermolecular, bool intramolecular, double? max_dist = null)
        {
            if (atoms == null)
            {
                throw  new ArgumentNullException(nameof(atoms));
            }

            var contact_element_types = new char[] {'C', 'N', 'O'};

            if (pdb_id != null) atoms = atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.OrdinalIgnoreCase)).ToList();

            if (intermolecular) atoms.ForEach(a => a.contact_map_intermolecular = new List<(atom atom, double distance)>());
            if (intramolecular) atoms.ForEach(a => a.contact_map_intramolecular = new List<(atom atom, double distance)>());

            if (contact_element_types != null && contact_element_types.Length > 0) atoms = atoms.Where(a => contact_element_types.Contains(a.element)).ToList();

            const bool skip_same_amino_acid_contacts = true;

            var pairs = atoms.SelectMany((a, x) => atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).Where(b => (intramolecular && b.atom1.chain_id == b.atom2.chain_id && (!skip_same_amino_acid_contacts || b.atom1.residue_index != b.atom2.residue_index)) || (intermolecular && b.atom1.chain_id != b.atom2.chain_id)).ToList()).ToList();

            var pair_distances = pairs.Select(a => (atom1: a.atom1, atom2: a.atom2, distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z))).ToList();

            for (var index = 0; index < pair_distances.Count; index++)
            {
                var a = pair_distances[index];

                if (max_dist != null && max_dist > 0 && a.distance > max_dist) continue;


                if (a.atom1.chain_id == a.atom2.chain_id)
                {
                    a.atom1.contact_map_intramolecular.Add((a.atom2, a.distance));
                    a.atom2.contact_map_intramolecular.Add((a.atom1, a.distance));
                }
                else //if (a.atom1.chain_id != a.atom2.chain_id)
                {
                    a.atom1.contact_map_intermolecular.Add((a.atom2, a.distance));
                    a.atom2.contact_map_intermolecular.Add((a.atom1, a.distance));
                }
            }
        }
        */

        //internal static List<(atom atom1, atom atom2, double distance)> get_master_to_master_contacts(string pdb_id, List<atom> atoms, string atom_type = null)
        //{
        //    if (atoms == null)
        //    {
        //        throw new ArgumentNullException(nameof(atoms));
        //    }

        //    var master_atoms = select_amino_acid_master_atoms(pdb_id, atoms, atom_type);

        //    var pairs = master_atoms.SelectMany((a, x) => master_atoms.Where((b, y) => x < y).Select((b, y) => (atom1: a, atom2: b)).ToList()).ToList();

        //    var pair_distances = pairs.Select(a => (atom1: a.atom1, atom2: a.atom2, distance: Distance3D(a.atom1.X, a.atom1.Y, a.atom1.Z, a.atom2.X, a.atom2.Y, a.atom2.Z))).ToList();

        //    return pair_distances;

        //}
    }
}