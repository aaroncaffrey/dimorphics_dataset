using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using dimorphics_dataset;
using RDotNet;

namespace protr_server
{
    public static class r_protr
    {
        public static readonly object engine_lock = new object();

        //public static REngine engine = null;//init_r();

        private static uint _key = 1;
        private static readonly object _key_lock = new object();

        public static int id;

        public static string Key
        {
            get
            {
                lock (_key_lock)
                {
                    _key++;
                    return $"{nameof(r_protr)}_{id}_{_key}";
                }
            }
        }

        //public r_peptides()
        //{
        //    engine = init_r();
        //}

        private static bool need_init = true;

        public static REngine init_r()
        {

            StartupParameter rinit = new StartupParameter();

            rinit.Quiet = true;
            rinit.Interactive = false;

            //rinit.RHome = $@"C:\Program Files\R\R-3.4.4\";
            rinit.RHome = $@"C:\Program Files\R\R-3.6.2\";

            //REngine.SetEnvironmentVariables();

            //C:\Program Files\R\R-3.6.2
            //var engine1 = REngine.GetInstance(@"C:\Program Files\R\R-3.4.4\bin\x64\R.dll", true, rinit);
            var engine1 = REngine.GetInstance(@"C:\Program Files\R\R-3.6.2\bin\x64\R.dll", true, rinit);

            if (need_init)
            {
                need_init = false;

                var r_init_cmds = $@"
                    #install.packages(""devtools"")
                    #library(devtools)
                    #install_github(""https://github.com/nanxstats/protr"")
                    library(protr)
                ";

                r_init_cmds.Split(new char[] { '\r', '\n' }).Where(a => !string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith("#")).ToList().ForEach(a => engine1.Evaluate(a));
            }

            return engine1;
        }

        public static List<subsequence_classification_data.feature_info> calculate_r_protr_classification_data_template = null;


        public static object get_values_lock = new object();

        public static List<subsequence_classification_data.feature_info> get_values(string seq, int id) //List<(string name, double value)>
        {
            lock (get_values_lock)
            {
                r_protr.id = id;

                //#if DEBUG
                ////if (program.verbose_debug) Program.WriteLine($"{nameof(get_values)}(string seq);");
                //#endif

                var engine = init_r();
                {



                    if (string.IsNullOrWhiteSpace(seq))
                    {
                        if (calculate_r_protr_classification_data_template != null)
                        {
                            var template = calculate_r_protr_classification_data_template.Select(a => new subsequence_classification_data.feature_info(a) { source = "", feature_value = 0 }).ToList();

                            return template;
                        }
                        else
                        {
                            throw new Exception();
                        }
                    }

                    var ret_extractAPAAC = extractAPAAC(engine, seq);
                    var ret_extractBLOSUM = extractBLOSUM(engine, seq);
                    
                    var ret_extractCTDC = extractCTDC(engine, seq);
                    var ret_extractCTDD = extractCTDD(engine, seq);
                    var ret_extractCTDT = extractCTDT(engine, seq);
                    var ret_extractCTriad = extractCTriad(engine, seq);

                    var ret_extractCTriadClass = extractCTriadClass(engine, seq);

                    var ret_extractDC = extractDC(engine, seq);

                    var ret_extractDescScales = extractDescScales(engine, seq);

                    var ret_extractFAScales = extractFAScales(engine, seq);

                    var ret_extractGeary = extractGeary(engine, seq);

                    var ret_extractMDSScales = extractMDSScales(engine, seq);

                    var ret_extractMoran = extractMoran(engine, seq);

                    var ret_extractMoreauBroto = extractMoreauBroto(engine, seq);
                    





                    
                    var features = new List<subsequence_classification_data.feature_info>();
                    var source = "";
                    var alpha_name = "Overall";

                    


                    if (calculate_r_protr_classification_data_template == null)
                    {
                        var template = features.Select(a => new subsequence_classification_data.feature_info(a) { source = "", feature_value = 0 }).ToList();
                        calculate_r_protr_classification_data_template = template;
                    }

                    return features;
                }
            }
        }


        public static List<(string name, int lambda, double value)> extractAPAAC(REngine engine, string x, int lambda_first=1, int lambda_last=2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(lambda_first), lambda_first.ToString()),
                (nameof(lambda_last), lambda_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractAPAAC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lambda, double value)>();

                for (var lambda = lambda_first; lambda <= lambda_last; lambda++)
                {
                    var f = nameof(extractAPAAC);
                    var k = Key;
                    var v = $"{f}_v{k}";
                    var ai = ""; //"[[1]]";

                    var evaluate = engine.Evaluate($"{v} <- {f}({nameof(x)} = \"{x}\", {nameof(lambda)} = {lambda})");
                    var names = engine.Evaluate($"names({v}{ai})").AsCharacter();
                    //var dimnames = engine.Evaluate($"dimnames({v}{ai})");
                    //var rownames = engine.Evaluate($"rownames({v}{ai})");
                    //var colnames = engine.Evaluate($"colnames({v}{ai})");
                    var values = engine.Evaluate($"{v}{ai}").AsNumeric();
                    var rm = engine.Evaluate($"rm({v})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) => (name: a, lambda: lambda, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                return list;
            }
        }



        public static List<(string submat, int k, int lag, string name, double value)> extractBLOSUM(REngine engine, string x, int lag_first =1, int lag_last=2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(lag_first), lag_first.ToString()),
                (nameof(lag_last), lag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractBLOSUM)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var submats = new string[]
                {
                    "AABLOSUM45", "AABLOSUM50", "AABLOSUM62", "AABLOSUM80", "AABLOSUM100", "AAPAM30", "AAPAM40",
                    "AAPAM70", "AAPAM120", "AAPAM250"
                };

                var k_first = 5;
                var k_last = 5;

                var list = new List<(string submat, int k, int lag, string name, double value)>();
                foreach (var submat in submats)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        for (var k = k_first; k <= k_last; k++)
                        {
                            var fn_name = nameof(extractBLOSUM);
                            var vr_name = $"{fn_name}_v{Key}";
                            var ar_ix = ""; //"[[1]]";

                            var evaluate =
                                engine.Evaluate(
                                    $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(submat)} = ""{submat}"", {nameof(k)} = {k}, {nameof(lag)} = {lag})");
                            var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                            //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                            //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                            //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                            var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                            var rm = engine.Evaluate($@"rm({vr_name})");

                            if (names.Length != values.Length) throw new Exception();

                            list.AddRange(names.Select((a, i) => (submat: submat, k: k, lag: lag, name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                        }
                    }
                }

                return list;
            }
        }


        public static List<(string name, double value)> extractCTDC(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractCTDC);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = ""; //"[[1]]";

                var evaluate = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }


        public static List<(string name, double value)> extractCTDD(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDD)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractCTDD);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = ""; //"[[1]]";

                var evaluate = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }



        public static List<(string name, double value)> extractCTDT(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDT)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractCTDT);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = ""; //"[[1]]";

                var evaluate = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }



        public static List<(string name, double value)> extractCTriad(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTriad)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractCTriad);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = ""; 

                var evaluate = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }
        

        public static List<(string name, double value)> extractCTriadClass(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTriadClass)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractCTriadClass);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = "";

                var aaclass = $"{fn_name}_v{Key}";

                var evaluate1 = engine.Evaluate($@"{aaclass} <- list( c(""G"", ""A"", ""S"", ""T"", ""P"", ""D"", ""C""), c(""N"", ""V"", ""E"", ""Q"", ""I"", ""L""), c(""M"", ""H"", ""K"", ""F"", ""R"", ""Y"", ""W"") )");
                var evaluate2 = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(aaclass)} = {aaclass})");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm1 = engine.Evaluate($@"rm({aaclass})");
                var rm2 = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }


        public static List<(string name, double value)> extractDC(REngine engine, string x)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractDC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractDC);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = ""; //"[[1]]";

                var evaluate = engine.Evaluate($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                return list;
            }
        }



        public static List<(string propmat, int pc, int lag, string name, List<(string row_name, string col_name, double pca_value)> pca, double value)> extractDescScales(REngine engine, string x, int lag_first=1, int lag_last=2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(lag_first), lag_first.ToString()),
                (nameof(lag_last), lag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractDescScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var propmats = new string[]
                {
                    "AAMOE2D", "AAMOE3D", "AACPSA", "AADescAll", "AA2DACOR", "AA3DMoRSE", "AAACF", "AABurden", "AAConn",
                    "AAConst", "AAEdgeAdj", "AAEigIdx", "AAFGC", "AAGeom", "AAGETAWAY", "AAInfo", "AAMolProp",
                    "AARandic", "AARDF", "AATopo", "AATopoChg", "AAWalk", "AAWHIM"
                };

                var index = "NULL";
                var silent = "FALSE";
                var scale = "TRUE";

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string propmat, int pc, int lag, string name, List<(string row_name, string col_name, double pca_value)> pca, double value)>();

                foreach (var propmat in propmats)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        for (var lag = 1; lag <= 2; lag++)
                        {
                            

                            var fn_name = nameof(extractDescScales);
                            var vr_name = $"{fn_name}_v{Key}";
                            var ar_ix = ""; //"[[1]]";

                            var vr_name_output = $"{fn_name}_v{Key}";

                            var cmd = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = ""{propmat}"", {nameof(index)} = {index}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";

                            var evaluate = engine.Evaluate(cmd);
                            var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                            //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                            //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                            //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                            var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                            var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacterMatrix();
                            var rm1 = engine.Evaluate($@"rm({vr_name})");
                            var rm2 = engine.Evaluate($@"rm({vr_name_output})");

                            var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                            var pca_col_labels = values_output[1, 0].Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);
                            for (var i = 2; i <= 4; i++) 
                            {
                                var pca_val = values_output[i,0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                                
                                var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                                var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                                pca_list.AddRange(pca_row_values.Select((a,j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                            } 

                            if (names.Length != values.Length) throw new Exception();

                            list.AddRange(names.Select((a, i) => (propmat: propmat, pc: pc, lag: lag, name: a, pca: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                        }
                    }
                }

                return list;
            }
        }



        public static List<(int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, string name, double value)> extractFAScales(REngine engine, string x, int lag_first=1, int lag_last=2, int factors_first=5, int factors_last=5)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(lag_first), lag_first.ToString()),
                (nameof(lag_last), lag_last.ToString()),
                (nameof(factors_first), factors_first.ToString()),
                (nameof(factors_last), factors_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractFAScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var silent = "FALSE";
                var scale = "TRUE";

                var list = new List<(int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, string name, double value)>();

                for (var factors = factors_first; factors <= factors_last; factors++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        var fn_name = nameof(extractFAScales);
                        var vr_name = $"{fn_name}_v{Key}";
                        var ar_ix = "";

                        var vr_name_tprops = $"{fn_name}_v{Key}";
                        var vr_name_output = $"{fn_name}_v{Key}";

                        var cmd1 = $@"{vr_name_tprops} <- AATopo[, c(37:41, 43:47)]"; // select a set of topological descriptors
                        var cmd2 = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}(x = ""{x}"", propmat = {vr_name_tprops}, {nameof(factors)} = {factors}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacter();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_tprops})");
                        var rm3 = engine.Evaluate($@"rm({vr_name_output})");

                        if (names.Length != values.Length) throw new Exception();

                        var uniqueness_labels = values_output[6].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var uniqueness_values = values_output[7].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                        // note: loadings_rows cannot be parsed correctly to be useful
                        var loadings_col_labels = values_output[10].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var loadings_rows = values_output.Skip(11).Take(uniqueness_labels.Length).Select(a => a.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)).ToList();

                        var factors_col_labels = values_output[22].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var factors_rows = values_output.Skip(23).Take(3).Select(a=>a.Split(new char[] {' '}, StringSplitOptions.RemoveEmptyEntries)).ToList();
                        var factors_row_labels = factors_rows.Select(a => string.Join(" ", a.Take(a.Length - factors_col_labels.Length).ToList())).ToList();
                        var factors_row_values = factors_rows.Select(a => a.Skip(a.Length - factors_col_labels.Length).Select(b=>double.Parse(b,NumberStyles.Float, CultureInfo.InvariantCulture)).ToList()).ToList();

                        var factors_list = new List<(string row_name, string col_name, double value)>();
                        for (var r = 0; r < factors_row_labels.Count; r++)
                        {
                            for (var c = 0; c < factors_col_labels.Length; c++)
                            {
                                factors_list.Add((factors_row_labels[r], factors_col_labels[c], factors_row_values[r][c]));
                            }
                        }
                        //var values_output2 = values_output[0,]

                        var chi_sq = double.Parse(values_output[28].Split()[5], NumberStyles.Float, CultureInfo.InvariantCulture);
                        var p_value = double.Parse(values_output[29].Split(new char[] { ' '}, StringSplitOptions.RemoveEmptyEntries).Last(),NumberStyles.Float,CultureInfo.InvariantCulture);

                        list.AddRange(names.Select((a, i) => (factors: factors, lag: lag, factors_list: factors_list, chi_sq: chi_sq, p_value: p_value, name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                return list;
            }
        }



        public static List<(string name, int nlag, double value)> extractGeary(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(nlag_first), nlag_first.ToString()),
                (nameof(nlag_last), nlag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractGeary)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    var fn_name = nameof(extractGeary);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = "";

                    var props = $@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"", ""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")";
                    var customprops = "NULL";

                    var cmd = $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) =>
                        (name: a, nlag: nlag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                return list;
            }
        }



        public static List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)> extractMDSScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(lag_first), lag_first.ToString()),
                (nameof(lag_last), lag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractMDSScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)>();

                var fn_name = nameof(extractMDSScales);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = "";

                var vr_name_tprops = $"{fn_name}_v{Key}";
                var vr_name_output = $"{fn_name}_v{Key}";

                var scale = "TRUE";
                var silent = "FALSE";

                var k_first = 5;
                var k_last = 5;

                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var k = k_first; k <= k_last; k++)
                    {
                        var cmd1 = $@"{vr_name_tprops} <- AATopo[, c(37:41, 43:47)]";
                        var cmd2 = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", propmat = {vr_name_tprops}, {nameof(k)} = {k}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";
                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacter();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_tprops})");
                        var rm3 = engine.Evaluate($@"rm({vr_name_output})");

                        var scaling_eigenvalues = values_output.Skip(1).SelectMany(a=> a.Split(new []{' '}, StringSplitOptions.RemoveEmptyEntries).Skip(1).Select(b=> double.Parse(b, NumberStyles.Float, CultureInfo.InvariantCulture)).ToArray()).ToArray();

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, k: k, lag: lag, scaling_eigenvalues: scaling_eigenvalues, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                return list;
            }
        }



        public static List<(string name, int nlag, double value)> extractMoran(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(nlag_first), nlag_first.ToString()),
                (nameof(nlag_last), nlag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractMoran)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    var fn_name = nameof(extractMoran);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = ""; //"[[1]]";

                    var props =
                        $@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"",""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")";
                    var customprops = "NULL";

                    var cmd =
                        $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) =>
                        (name: a, nlag: nlag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                return list;
            }
        }




        public static List<(string name, int nlag, double value)> extractMoreauBroto(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine.ToString()),
                (nameof(x), x.ToString()),
                (nameof(nlag_first), nlag_first.ToString()),
                (nameof(nlag_last), nlag_last.ToString()),
            };
            Console.WriteLine($@"{nameof(r_protr)}.{nameof(extractMoreauBroto)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();
                
                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    var fn_name = nameof(extractMoran);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = ""; //"[[1]]";

                
                    var props =
                        $@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"",""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")";
                    var customprops = "NULL";

                    var cmd =
                        $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) =>
                        (name: a, nlag: nlag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                return list;
            }



        }

    }
}
