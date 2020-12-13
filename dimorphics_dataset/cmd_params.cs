using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;

namespace dimorphics_dataset
{

    internal class cmd_params
    {
        //internal string[] area;
        internal bool use_dssp3;
        internal int class_id;
        internal string class_name;
        internal int min_sequence_length;
        internal int max_features;
        internal string output_folder;
        internal int? first_index;
        internal int? last_index;
        internal bool? verbose;
        //internal string tag;
        internal bool? use_children;

        internal cmd_params()
        {

        }

        internal cmd_params(cmd_params p)
        {
            //area = p.area;
            use_dssp3 = p.use_dssp3;
            class_id = p.class_id;
            class_name = /*program.string_debug*/(p.class_name);
            min_sequence_length = p.min_sequence_length;
            max_features = p.max_features;
            output_folder = /*program.string_debug*/(p.output_folder);
            first_index = p.first_index;
            last_index = p.last_index;
            verbose = p.verbose;
            use_children = p.use_children;
            //tag = p.tag;
        }

        internal string get_output_file_tag()
        {
            var x = new List<string>();

            var feature_types_1d_interface_kvl = feature_types_1d_interface.key_value_list();
            if (feature_types_1d_interface_kvl.All(a => a.value))
            {
                x.Add(area1i);
            }
            else if (feature_types_1d_interface_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_1d_interface_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area1i}.{a.key}")).ToArray());
            }

            var feature_types_1d_neighbourhood_kvl = feature_types_1d_neighbourhood.key_value_list();
            if (feature_types_1d_neighbourhood_kvl.All(a => a.value))
            {
                x.Add(area1n);
            }
            else if (feature_types_1d_neighbourhood_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_1d_neighbourhood_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area1n}.{a.key}")).ToArray());
            }

            var feature_types_1d_chain_kvl = feature_types_1d_chain.key_value_list();
            if (feature_types_1d_chain_kvl.All(a => a.value))
            {
                x.Add(area1p);
            }
            else if (feature_types_1d_chain_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_1d_chain_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area1p}.{a.key}")).ToArray());
            }



            var feature_types_2d_interface_kvl = feature_types_2d_interface.key_value_list();
            if (feature_types_2d_interface_kvl.All(a => a.value))
            {
                x.Add(area2i);
            }
            else if (feature_types_2d_interface_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_2d_interface_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area2i}.{a.key}")).ToArray());
            }

            var feature_types_2d_neighbourhood_kvl = feature_types_2d_neighbourhood.key_value_list();
            if (feature_types_2d_neighbourhood_kvl.All(a => a.value))
            {
                x.Add(area2n);
            }
            else if (feature_types_2d_neighbourhood_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_2d_neighbourhood_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area2n}.{a.key}")).ToArray());
            }

            var feature_types_2d_chain_kvl = feature_types_2d_chain.key_value_list();
            if (feature_types_2d_chain_kvl.All(a => a.value))
            {
                x.Add(area2p);
            }
            else if (feature_types_2d_chain_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_2d_chain_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area2p}.{a.key}")).ToArray());
            }



            var feature_types_3d_interface_kvl = feature_types_3d_interface.key_value_list();
            if (feature_types_3d_interface_kvl.All(a => a.value))
            {
                x.Add(area3i);
            }
            else if (feature_types_3d_interface_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_3d_interface_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area3i}.{a.key}")).ToArray());
            }

            var feature_types_3d_neighbourhood_kvl = feature_types_3d_neighbourhood.key_value_list();
            if (feature_types_3d_neighbourhood_kvl.All(a => a.value))
            {
                x.Add(area3n);
            }
            else if (feature_types_3d_neighbourhood_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_3d_neighbourhood_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area3n}.{a.key}")).ToArray());
            }

            var feature_types_3d_chain_kvl = feature_types_3d_chain.key_value_list();
            if (feature_types_3d_chain_kvl.All(a => a.value))
            {
                x.Add(area3p);
            }
            else if (feature_types_3d_chain_kvl.Any(a => a.value))
            {
                x.AddRange(feature_types_3d_chain_kvl.Where(a => a.value).Select(a => /*program.string_debug*/($"{area3p}.{a.key}")).ToArray());
            }

            // f_([1i]_[1n.aaindex]_[1n.dfdsf]_[1p])_(+1)_(dimorphic_coil).csv

            return string.Join(/*program.string_debug*/($@"_"), x.Select(a => /*program.string_debug*/($"[{a}]")).ToArray());
        }


        internal static List<(string key, string value)> get_params(string[] args)
        {
            var x = new List<(string key, string value)>();

            if (args == null || args.Length == 0)
            {
                return x;
            }

            args = args.SelectMany(a => a.Split(new char[] { '=' }, StringSplitOptions.RemoveEmptyEntries)).Where(a => !string.IsNullOrWhiteSpace(a)).ToArray();

            var name = "";
            var value = "";

            for (var i = 0; i < args.Length; i++)
            {
                var arg = args[i];

                var starts_dash = arg[0] == '-';
                var is_num = double.TryParse(arg, NumberStyles.Float, NumberFormatInfo.InvariantInfo, out var val_out);

                var is_name = starts_dash && (!is_num || string.IsNullOrEmpty(name)) && arg.Length > 1;
                var is_value = !is_name;
                var is_final_index = i == args.Length - 1;

                if (is_name)
                {
                    if (!string.IsNullOrWhiteSpace(name) || !string.IsNullOrWhiteSpace(value))
                    {
                        x.Add((name, value));
                    };

                    name = arg[1..];
                    value = "";
                }

                if (is_value) value += (value.Length > 0 ? " " : "") + arg;

                if (is_final_index)
                {
                    x.Add((name, value));
                }
            }

            return x;
        }

        internal const string area1i = @"1i";
        internal const string area1n = @"1n";
        internal const string area1p = @"1p";

        internal const string area2i = @"2i";
        internal const string area2n = @"2n";
        internal const string area2p = @"2p";

        internal const string area3i = @"3i";
        internal const string area3n = @"3n";
        internal const string area3p = @"3p";

        internal static readonly string[] defined_areas = { area1i, area1n, area1p, area2i, area2n, area2p, area3i, area3n, area3p };

        internal feature_types_1d feature_types_1d_interface = new feature_types_1d(enum_protein_data_source.interface_1d, false);
        internal feature_types_1d feature_types_1d_neighbourhood = new feature_types_1d(enum_protein_data_source.neighbourhood_1d, false);
        internal feature_types_1d feature_types_1d_chain = new feature_types_1d(enum_protein_data_source.chain_1d, false);

        internal feature_types_2d feature_types_2d_interface = new feature_types_2d(enum_protein_data_source.interface_2d, false);
        internal feature_types_2d feature_types_2d_neighbourhood = new feature_types_2d(enum_protein_data_source.neighbourhood_2d, false);
        internal feature_types_2d feature_types_2d_chain = new feature_types_2d(enum_protein_data_source.chain_2d, false);

        internal feature_types_3d feature_types_3d_interface = new feature_types_3d(enum_protein_data_source.interface_3d, false);
        internal feature_types_3d feature_types_3d_neighbourhood = new feature_types_3d(enum_protein_data_source.neighbourhood_3d, false);
        internal feature_types_3d feature_types_3d_chain = new feature_types_3d(enum_protein_data_source.chain_3d, false);



        internal cmd_params(string[] args)
        {
            const string module_name = nameof(cmd_params);
            const string method_name = nameof(cmd_params);

            const string standard_coil = @"standard_coil";
            const string dimorphic_coil = @"dimorphic_coil";

            var param_list = get_params(args);
            io_proxy.WriteLine(string.Join(/*program.string_debug*/($@", "), param_list), nameof(program), nameof(cmd_params));

            var known_params = new List<string>()
            {
                //nameof(area),
                nameof(use_dssp3),
                nameof(class_id),
                nameof(class_name),
                nameof(min_sequence_length),
                nameof(max_features),
                nameof(output_folder),
                nameof(first_index),
                nameof(last_index),
                nameof(verbose),
                //nameof(tag),
                nameof(use_children),
            };

            known_params.AddRange(defined_areas);

            known_params.AddRange(feature_types_1d.keys.Select(a => /*program.string_debug*/($"{area1i}.{a}")).ToArray());
            known_params.AddRange(feature_types_1d.keys.Select(a => /*program.string_debug*/($"{area1n}.{a}")).ToArray());
            known_params.AddRange(feature_types_1d.keys.Select(a => /*program.string_debug*/($"{area1p}.{a}")).ToArray());
                                                                    
            known_params.AddRange(feature_types_2d.keys.Select(a => /*program.string_debug*/($"{area2i}.{a}")).ToArray());
            known_params.AddRange(feature_types_2d.keys.Select(a => /*program.string_debug*/($"{area2n}.{a}")).ToArray());
            known_params.AddRange(feature_types_2d.keys.Select(a => /*program.string_debug*/($"{area2p}.{a}")).ToArray());
                                                                    
            known_params.AddRange(feature_types_3d.keys.Select(a => /*program.string_debug*/($"{area3i}.{a}")).ToArray());
            known_params.AddRange(feature_types_3d.keys.Select(a => /*program.string_debug*/($"{area3n}.{a}")).ToArray());
            known_params.AddRange(feature_types_3d.keys.Select(a => /*program.string_debug*/($"{area3p}.{a}")).ToArray());

            var unknown_params = param_list.Where(a => known_params.All(b => !string.Equals(a.key, b, StringComparison.OrdinalIgnoreCase))).Select(a => a.key).ToArray();

            if (unknown_params != null && unknown_params.Length > 0)
            {
                throw new ArgumentOutOfRangeException(nameof(args), /*program.string_debug*/($@"Arguments not recognised: {string.Join(/*program.string_debug*/($@", "), unknown_params)}"));
            }


            if (param_list.Count == 0)
            {
                var classes = new[] { (class_id: +1, class_name: dimorphic_coil), (class_id: -1, class_name: standard_coil) };
                var output_folder = /*program.string_debug*/(@"e:\dataset\");
                program.verbose = true;
                var exe = Path.GetFileName(Process.GetCurrentProcess().MainModule.FileName);
                io_proxy.WriteLine(/*program.string_debug*/($@""));
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(dimorphics_dataset)} usage:"));
                io_proxy.WriteLine(/*program.string_debug*/($@""));

                foreach (var defined_area in defined_areas)
                {
                    io_proxy.WriteLine(/*program.string_debug*/($@"{defined_area[0]}d {(defined_area[1] == 'i' || defined_area[1] == 'I' ? /*program.string_debug*/($@"interface subsequence") : /*program.string_debug*/($@""))}{(defined_area[1] == 'n' || defined_area[1] == 'N' ? /*program.string_debug*/($@"neighbourhood") : /*program.string_debug*/($@""))}{(defined_area[1] == 'p' || defined_area[1] == 'P' ? /*program.string_debug*/($@"protein chain") : /*program.string_debug*/($@""))} area:"));

                    foreach (var b in classes)
                    {
                        var stdout_file = Path.Combine(output_folder, /*program.string_debug*/($@"stdout_{defined_area}_({(b.class_id > 0 ? /*program.string_debug*/($@"+") : /*program.string_debug*/($@""))}{b.class_id})_({b.class_name}).txt"));
                        var stderr_file = Path.Combine(output_folder, /*program.string_debug*/($@"stderr_{defined_area}_({(b.class_id > 0 ? /*program.string_debug*/($@"+") : /*program.string_debug*/($@""))}{b.class_id})_({b.class_name}).txt"));

                        string[] keys;
                        if (new string[] { area1i, area1n, area1p }.Any(a => string.Equals(a, defined_area, StringComparison.OrdinalIgnoreCase))) keys = feature_types_1d.keys.Select(a => /*program.string_debug*/($"{defined_area}.{a}")).ToArray();
                        else if (new string[] { area2i, area2n, area2p }.Any(a => string.Equals(a, defined_area, StringComparison.OrdinalIgnoreCase))) keys = feature_types_2d.keys.Select(a => /*program.string_debug*/($"{defined_area}.{a}")).ToArray();
                        else if (new string[] { area3i, area3n, area3p }.Any(a => string.Equals(a, defined_area, StringComparison.OrdinalIgnoreCase))) keys = feature_types_3d.keys.Select(a => /*program.string_debug*/($"{defined_area}.{a}")).ToArray();
                        else throw new Exception();

                        io_proxy.WriteLine(/*program.string_debug*/($@"{exe} -{defined_area}=true {string.Join(/*program.string_debug*/($@" "), keys.Select(a => /*program.string_debug*/($"-{a}=true")).ToArray())} -{nameof(use_dssp3)}=true -{nameof(class_id)}={(b.class_id > 0 ? /*program.string_debug*/($@"+") : /*program.string_debug*/($@""))}{b.class_id} -{nameof(class_name)}={b.class_name} -{nameof(min_sequence_length)}=3 -{nameof(max_features)}=100 -{nameof(output_folder)}={output_folder} -{nameof(use_children)}=true -{nameof(verbose)}=true [-{nameof(first_index)}=[i] -{nameof(last_index)}=[j]] 1> {stdout_file} 2> {stderr_file}"));

                    }

                    io_proxy.WriteLine(/*program.string_debug*/($@""));
                }

                Environment.Exit(0);
                //return default;
                //throw new Exception(/*program.string_debug*/($@"No {nameof(args)} specified.");
            }

            if (param_list.Any(a => string.Equals(a.key, area1i, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1i, StringComparison.OrdinalIgnoreCase));
                feature_types_1d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area1n, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1n, StringComparison.OrdinalIgnoreCase));
                feature_types_1d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area1p, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area1p, StringComparison.OrdinalIgnoreCase));
                feature_types_1d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2i, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2i, StringComparison.OrdinalIgnoreCase));
                feature_types_2d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2n, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2n, StringComparison.OrdinalIgnoreCase));
                feature_types_2d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area2p, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area2p, StringComparison.OrdinalIgnoreCase));
                feature_types_2d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3i, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3i, StringComparison.OrdinalIgnoreCase));
                feature_types_3d_interface.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3n, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3n, StringComparison.OrdinalIgnoreCase));
                feature_types_3d_neighbourhood.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            if (param_list.Any(a => string.Equals(a.key, area3p, StringComparison.OrdinalIgnoreCase)))
            {
                var v = param_list.First(a => string.Equals(a.key, area3p, StringComparison.OrdinalIgnoreCase));
                feature_types_3d_chain.set_enable(string.IsNullOrEmpty(v.value) || string.Equals(v.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(v.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase));
            }

            ///*program.string_debug*/()

            feature_types_1d_interface
              .set_enable(param_list
                  .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area1i}."), StringComparison.OrdinalIgnoreCase) && feature_types_1d.keys.Contains(a.key[/*program.string_debug*/($"{area1i}.").Length..])
                  )
                  .Select(a =>
                      (
                          enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                          name: a.key[/*program.string_debug*/($"{area1i}.").Length..]
                      )
                  )
                  .ToArray()
              );

            feature_types_1d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area1n}."), StringComparison.OrdinalIgnoreCase) && feature_types_1d.keys.Contains(a.key[$"{area1n}.".Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area1n}.").Length..]
                        )
                    )
                    .ToArray()
                );


            feature_types_1d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area1p}."), StringComparison.OrdinalIgnoreCase) && feature_types_1d.keys.Contains(a.key[/*program.string_debug*/($"{area1p}.").Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area1p}.").Length..]
                        )
                    )
                    .ToArray()
                );

            feature_types_2d_interface
                     .set_enable(param_list
                         .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area2i}."), StringComparison.OrdinalIgnoreCase) && feature_types_2d.keys.Contains(a.key[/*program.string_debug*/($"{area2i}.").Length..])
                         )
                         .Select(a =>
                             (
                                 enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, "true", StringComparison.OrdinalIgnoreCase)),
                                 name: a.key[/*program.string_debug*/($"{area2i}.").Length..]
                             )
                         )
                         .ToArray()
                     );

            feature_types_2d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area2n}."), StringComparison.OrdinalIgnoreCase) && feature_types_2d.keys.Contains(a.key[$"{area2n}.".Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, "1", StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, "true", StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area2n}.").Length..]
                        )
                    )
                    .ToArray()
                );


            feature_types_2d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area2p}."), StringComparison.OrdinalIgnoreCase) && feature_types_2d.keys.Contains(a.key[/*program.string_debug*/($"{area2p}.").Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area2p}.").Length..]
                        )
                    )
                    .ToArray()
                );


            feature_types_3d_interface
                     .set_enable(param_list
                         .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area3i}."), StringComparison.OrdinalIgnoreCase) && feature_types_3d.keys.Contains(a.key[/*program.string_debug*/($"{area3i}.").Length..])
                         )
                         .Select(a =>
                             (
                                 enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                                 name: a.key[/*program.string_debug*/($"{area3i}.").Length..]
                             )
                         )
                         .ToArray()
                     );

            feature_types_3d_neighbourhood
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area3n}."), StringComparison.OrdinalIgnoreCase) && feature_types_3d.keys.Contains(a.key[/*program.string_debug*/($"{area3n}.").Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area3n}.").Length..]
                        )
                    )
                    .ToArray()
                );


            feature_types_3d_chain
                .set_enable(param_list
                    .Where(a => a.key.StartsWith(/*program.string_debug*/($"{area3p}."), StringComparison.OrdinalIgnoreCase) && feature_types_3d.keys.Contains(a.key[/*program.string_debug*/($"{area3p}.").Length..])
                    )
                    .Select(a =>
                        (
                            enable: (string.IsNullOrEmpty(a.value) || string.Equals(a.value, /*program.string_debug*/("1"), StringComparison.OrdinalIgnoreCase) || string.Equals(a.value, /*program.string_debug*/("true"), StringComparison.OrdinalIgnoreCase)),
                            name: a.key[/*program.string_debug*/($"{area3p}.").Length..]
                        )
                    )
                    .ToArray()
                );





            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(use_dssp3), StringComparison.OrdinalIgnoreCase)).value))
            {
                use_dssp3 = bool.Parse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(use_dssp3), StringComparison.OrdinalIgnoreCase)).value);
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(use_dssp3)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(class_id), StringComparison.OrdinalIgnoreCase)).value))
            {
                class_id = int.Parse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(class_id), StringComparison.OrdinalIgnoreCase)).value, NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(class_id)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(class_name), StringComparison.OrdinalIgnoreCase)).value))
            {
                class_name = /*program.string_debug*/(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(class_name), StringComparison.OrdinalIgnoreCase)).value?.ToLowerInvariant());
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(class_name)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(min_sequence_length), StringComparison.OrdinalIgnoreCase)).value))
            {
                min_sequence_length = int.Parse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(min_sequence_length), StringComparison.OrdinalIgnoreCase)).value, NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(min_sequence_length)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(max_features), StringComparison.OrdinalIgnoreCase)).value))
            {
                max_features = int.Parse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(max_features), StringComparison.OrdinalIgnoreCase)).value, NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(max_features)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(output_folder), StringComparison.OrdinalIgnoreCase)).value))
            {
                output_folder = /*program.string_debug*/(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(output_folder), StringComparison.OrdinalIgnoreCase)).value);
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(output_folder)}"));
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(first_index), StringComparison.OrdinalIgnoreCase)).value))
            {
                first_index = int.TryParse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(first_index), StringComparison.OrdinalIgnoreCase)).value, NumberStyles.Integer, NumberFormatInfo.InvariantInfo, out var tp_first_index) ? tp_first_index : (int?)null;
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(last_index), StringComparison.OrdinalIgnoreCase)).value))
            {
                last_index = int.TryParse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(last_index), StringComparison.OrdinalIgnoreCase)).value, NumberStyles.Integer, NumberFormatInfo.InvariantInfo, out var tp_last_index) ? tp_last_index : (int?)null;
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(verbose), StringComparison.OrdinalIgnoreCase)).value))
            {
                verbose = bool.TryParse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(verbose), StringComparison.OrdinalIgnoreCase)).value, out var tp_verbose) ? tp_verbose : (bool?)null;
            }

            if (!string.IsNullOrWhiteSpace(param_list.FirstOrDefault(a => string.Equals(a.key, nameof(use_children), StringComparison.OrdinalIgnoreCase)).value))
            {
                use_children = bool.TryParse(param_list.FirstOrDefault(z => string.Equals(z.key, nameof(use_children), StringComparison.OrdinalIgnoreCase)).value, out var tp_use_children) ? tp_use_children : (bool?)null;
            }
            else
            {
                throw new ArgumentNullException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: missing argument {nameof(use_children)}"));
            }


            verbose ??= (first_index == null && last_index == null);
            first_index ??= last_index;
            last_index ??= first_index;
            use_children ??= true;


            if (first_index == null && last_index == null)
            {
                io_proxy.WriteLine();

                //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(area)} = ""{string.Join(/*program.string_debug*/($@", "), area ?? Array.Empty<string>())}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(use_dssp3)} = ""{use_dssp3}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(class_id)} = ""{class_id}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(class_name)} = ""{class_name}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(min_sequence_length)} = ""{min_sequence_length}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(max_features)} = ""{max_features}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(output_folder)} = ""{output_folder}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(first_index)} = ""{first_index}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(last_index)} = ""{last_index}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(verbose)} = ""{verbose}"""), module_name, method_name);
                io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(use_children)} = ""{use_children}"""), module_name, method_name);

                io_proxy.WriteLine();

                //feature_types_1d_interface.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_1d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_1d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_1d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_1d_chain.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_1d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_2d_interface.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_2d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_2d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_2d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_2d_chain.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_2d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_3d_interface.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_3d_interface)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_3d_neighbourhood.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_3d_neighbourhood)}.{a.key} = ""{a.value}""", module_name, method_name));

                //io_proxy.WriteLine();

                //feature_types_3d_chain.key_value_list().ForEach(a => io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(feature_types_3d_chain)}.{a.key} = ""{a.value}""", module_name, method_name));
            }

            var feature_status = new List<(string key, bool value)>();

            feature_status.AddRange(feature_types_1d_interface.key_value_list().Select(a => (/*program.string_debug*/($"{area1i}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_1d_neighbourhood.key_value_list().Select(a => (/*program.string_debug*/($"{area1n}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_1d_chain.key_value_list().Select(a => (/*program.string_debug*/($"{area1p}.{a.key}"), a.value)).ToArray());

            feature_status.AddRange(feature_types_2d_interface.key_value_list().Select(a => (/*program.string_debug*/($"{area2i}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_2d_neighbourhood.key_value_list().Select(a => (/*program.string_debug*/($"{area2n}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_2d_chain.key_value_list().Select(a => (/*program.string_debug*/($"{area2p}.{a.key}"), a.value)).ToArray());

            feature_status.AddRange(feature_types_3d_interface.key_value_list().Select(a => (/*program.string_debug*/($"{area3i}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_3d_neighbourhood.key_value_list().Select(a => (/*program.string_debug*/($"{area3n}.{a.key}"), a.value)).ToArray());
            feature_status.AddRange(feature_types_3d_chain.key_value_list().Select(a => (/*program.string_debug*/($"{area3p}.{a.key}"), a.value)).ToArray());


            var feature_enabled = feature_status.Where(a => a.value).ToList();
            var num_enabled = feature_enabled.Count;
            io_proxy.WriteLine(/*program.string_debug*/($@"Features enabled for extraction [{num_enabled}]:"), module_name, method_name);
            feature_enabled.ForEach(a => io_proxy.WriteLine(a.key, module_name, method_name));
            io_proxy.WriteLine();

            var feature_disabled = feature_status.Where(a => !a.value).ToList();
            var num_disabled = feature_disabled.Count;
            io_proxy.WriteLine(/*program.string_debug*/($@"Features disabled for extraction [{num_disabled}]:"), module_name, method_name);
            feature_disabled.ForEach(a => io_proxy.WriteLine(a.key, module_name, method_name));
            io_proxy.WriteLine();
        

            if (num_enabled == 0)
            {
                throw new ArgumentOutOfRangeException(nameof(args), /*program.string_debug*/($@"{module_name}.{method_name}: No features are selected for extraction ({nameof(num_enabled)}: {num_enabled}, {nameof(num_disabled)}: {num_disabled})."));
            }
        }
    }

}
/*internal static feature_types feature_types_params(string[] area)
{
    const string module_name = nameof(dimorphics_dataset.feature_types);
    const string method_name = nameof(feature_types_params);

    if (area == null || area.Length == 0)
    {
        return null;
    }

    area = area.SelectMany(a => a.Split(new char[] { ' ', ';', ',', '-', '|' }, StringSplitOptions.RemoveEmptyEntries)).Where(a => !string.IsNullOrWhiteSpace(a)).ToArray();


    var do_1d_interface = area?.Any(a => string.Equals(a, cmd_params.area1i, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_1d_nh = area?.Any(a => string.Equals(a, cmd_params.area1n, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_1d_protein = area?.Any(a => string.Equals(a, cmd_params.area1p, StringComparison.OrdinalIgnoreCase)) ?? false;

    var do_2d_interface = area?.Any(a => string.Equals(a, cmd_params.area2i, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_2d_nh = area?.Any(a => string.Equals(a, cmd_params.area2n, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_2d_protein = area?.Any(a => string.Equals(a, cmd_params.area2p, StringComparison.OrdinalIgnoreCase)) ?? false;

    var do_3d_interface = area?.Any(a => string.Equals(a, cmd_params.area3i, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_3d_nh = area?.Any(a => string.Equals(a, cmd_params.area3n, StringComparison.OrdinalIgnoreCase)) ?? false;
    var do_3d_protein = area?.Any(a => string.Equals(a, cmd_params.area3p, StringComparison.OrdinalIgnoreCase)) ?? false;


    var feature_types = new feature_types()
    {
        // 1d
        feature_types_1d_interface = !do_1d_interface
            ? null
            : new feature_types_1d(enum_protein_data_source.interface_1d)
            {
                pse_aac = do_1d_interface,
                geometry = do_1d_interface,
                        //mpsa_classification_data_subsequence = do_1d_interface,
                        iupred2a = do_1d_interface,
                aaindex = do_1d_interface,
                sable = do_1d_interface,
                stackdppred = false, //must be false - 1d protein level only
                        blast_pssm = do_1d_interface,
                r_peptides = do_1d_interface,
                r_protr = do_1d_interface,
            },
        feature_types_neighbourhood_1d = !do_1d_nh
            ? null
            : new feature_types_1d(enum_protein_data_source.neighbourhood_1d)
            {
                pse_aac = do_1d_nh,
                geometry = do_1d_nh,
                        //mpsa_classification_data_subsequence = do_1d_nh,
                        iupred2a = do_1d_nh,
                aaindex = do_1d_nh,
                sable = do_1d_nh,
                stackdppred = false, //must be false - 1d protein level only
                        blast_pssm = do_1d_nh,
                r_peptides = do_1d_nh,
                r_protr = do_1d_nh,
            },
        feature_types_chain_1d = !do_1d_protein
            ? null
            : new feature_types_1d(enum_protein_data_source.chain_1d)
            {
                pse_aac = do_1d_protein,
                geometry = do_1d_protein,
                        //mpsa_classification_data_subsequence = do_1d_protein,
                        iupred2a = do_1d_protein,
                aaindex = do_1d_protein,
                sable = do_1d_protein,
                stackdppred = do_1d_protein,
                blast_pssm = do_1d_protein,
                r_peptides = do_1d_protein,
                r_protr = do_1d_protein,
            },

        // 2d
        feature_types_2d_interface = !do_2d_interface
            ? null
            : new feature_types_2d(enum_protein_data_source.interface_2d)
            {
                mpsa = do_2d_interface,
            },
        feature_types_neighbourhood_2d = !do_2d_nh
            ? null
            : new feature_types_2d(enum_protein_data_source.neighbourhood_2d)
            {

                mpsa = do_2d_nh,

            },
        feature_types_chain_2d = !do_2d_protein
            ? null
            : new feature_types_2d(enum_protein_data_source.chain_2d)
            {
                mpsa = do_2d_protein,

            },


        // 3d

        feature_types_3d_interface = !do_3d_interface
            ? null
            : new feature_types_3d(enum_protein_data_source.interface_3d)
            {
                sasa = do_3d_interface,
                intramolecular = do_3d_interface,

                foldx = do_3d_interface,
                tortuosity = do_3d_interface,
                ring = do_3d_interface,
                pse_ssc_dssp = do_3d_interface,
            },
        feature_types_neighbourhood_3d = !do_3d_nh
            ? null
            : new feature_types_3d(enum_protein_data_source.neighbourhood_3d)
            {
                sasa = do_3d_nh,
                intramolecular = do_3d_nh,

                foldx = do_3d_nh,
                tortuosity = do_3d_nh,
                ring = do_3d_nh,
                pse_ssc_dssp = do_3d_nh,
            },
        feature_types_chain_3d = !do_3d_protein
            ? null
            : new feature_types_3d(enum_protein_data_source.chain_3d)
            {
                sasa = do_3d_protein,
                intramolecular = do_3d_protein,

                foldx = do_3d_protein,
                tortuosity = do_3d_protein,
                ring = do_3d_protein,
                pse_ssc_dssp = do_3d_protein,
            }
    };

    // 1d                               
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_1d_interface)} = {feature_types.feature_types_1d_interface}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_neighbourhood_1d)} = {feature_types.feature_types_neighbourhood_1d}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_chain_1d)} = {feature_types.feature_types_chain_1d}"), module_name, method_name);
                                        
    //2d                                
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_2d_interface)} = {feature_types.feature_types_2d_interface}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_neighbourhood_2d)} = {feature_types.feature_types_neighbourhood_2d}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_chain_2d)} = {feature_types.feature_types_chain_2d}"), module_name, method_name);
                                        
    //3d                                
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_3d_interface)} = {feature_types.feature_types_3d_interface}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}"), module_name, method_name);
    io_proxy.WriteLine(/*program.string_debug* /($@"{nameof(feature_types.feature_types_chain_3d)} = {feature_types.feature_types_chain_3d}"), module_name, method_name);
    return feature_types;
}*/