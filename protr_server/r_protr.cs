using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using dimorphics_dataset;
using RDotNet;

namespace protr_server
{
    public static class r_protr
    {
        public static readonly object engine_lock = new object();
        private static uint _key = 1;
        private static readonly object _key_lock = new object();
        private static int _id;
        private static List<feature_info> _template_get_values = null;
        private static readonly object _get_values_lock = new object();

        public static string Key
        {
            get
            {
                lock (_key_lock)
                {
                    _key++;
                    return $"{nameof(r_protr)}_{_id}_{_key}";
                }
            }
        }

        private static bool need_init = true;

        public static REngine init_r()
        {
            StartupParameter rinit = new StartupParameter();
            rinit.Quiet = true;
            rinit.Interactive = false;
            rinit.RHome = $@"C:\Program Files\R\R-3.6.2\";
            var engine1 = REngine.GetInstance(Path.Combine(rinit.RHome, $@"bin\x64\R.dll"), true, rinit);

            if (need_init)
            {
                need_init = false;

                var r_init_cmds = $@"
                    #install.packages(""devtools"")
                    #library(devtools)
                    #install_github(""https://github.com/nanxstats/protr"")
                    library(protr)
                ";

                r_init_cmds.Split(new char[] { '\r', '\n' }).Where(a => !string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith("#", StringComparison.InvariantCulture)).ToList().ForEach(a => engine1.Evaluate(a));
            }

            return engine1;
        }

        public static List<feature_info> get_values(int id, string source, string alphabet_name, string sequence)
        {
            const int min_sequence_length = 3;

            lock (_get_values_lock)
            {
                r_protr._id = id;

                var engine = init_r();
                {
                    if (string.IsNullOrWhiteSpace(sequence) || sequence.Length < min_sequence_length)
                    {
                        if (_template_get_values == null)
                        {
                            _template_get_values = get_values(id, source, alphabet_name, "ALG");
                        }

                        if (_template_get_values != null)
                        {
                            var template = _template_get_values.Select(a => new feature_info(a) { alphabet = alphabet_name, source = source, feature_value = 0 }).ToList();

                            return template;
                        }

                        throw new Exception();
                    }

                    var ret_extractAPAAC = extractAPAAC(engine, sequence);
                    var ret_extractBLOSUM = extractBLOSUM(engine, sequence);
                    var ret_extractCTDC = extractCTDC(engine, sequence);
                    var ret_extractCTDD = extractCTDD(engine, sequence);
                    var ret_extractCTDT = extractCTDT(engine, sequence);
                    var ret_extractCTriad = extractCTriad(engine, sequence);
                    var ret_extractCTriadClass = extractCTriadClass(engine, sequence);
                    var ret_extractDC = extractDC(engine, sequence);
                    var ret_extractDescScales = extractDescScales(engine, sequence);
                    var ret_extractFAScales = extractFAScales(engine, sequence);
                    var ret_extractGeary = extractGeary(engine, sequence);
                    var ret_extractMDSScales = extractMDSScales(engine, sequence);
                    var ret_extractMoran = extractMoran(engine, sequence);
                    var ret_extractMoreauBroto = extractMoreauBroto(engine, sequence);
                    var ret_extractPAAC = extractPAAC(engine, sequence);
                    var ret_extractProtFP = extractProtFP(engine, sequence);
                    var ret_extractProtFPGap = extractProtFPGap(engine, sequence);
                    var ret_extractQSO = extractQSO(engine, sequence);
                    var ret_extractSOCN = extractSOCN(engine, sequence);
                    var ret_extractScales = extractScales(engine, sequence);
                    var ret_extractScalesGap = extractScalesGap(engine, sequence);
                    //var ret_extractTC =  extractTC(engine, sequence); // 800 features


                    var features = new List<feature_info>();

                    var ret_extractAPAAC_f = ret_extractAPAAC.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractAPAAC)}_{a.lambda}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractBLOSUM_f = ret_extractBLOSUM.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractBLOSUM)}_{a.submat}_{a.k}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractCTDC_f = ret_extractCTDC.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractCTDC)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractCTDD_f = ret_extractCTDD.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractCTDD)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractCTDT_f = ret_extractCTDT.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractCTDT)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractCTriad_f = ret_extractCTriad.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractCTriad)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractCTriadClass_f = ret_extractCTriadClass.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractCTriadClass)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractDC_f = ret_extractDC.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractDC)}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractDescScales_f = ret_extractDescScales.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractDescScales)}_{a.propmat}_{a.pc}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractFAScales_f = ret_extractFAScales.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractFAScales)}_{a.factors}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractGeary_f = ret_extractGeary.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractGeary)}_{a.nlag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractMDSScales_f = ret_extractMDSScales.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractMDSScales)}_{a.k}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractMoran_f = ret_extractMoran.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractMoran)}_{a.nlag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractMoreauBroto_f = ret_extractMoreauBroto.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractMoreauBroto)}_{a.nlag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractPAAC_f = ret_extractPAAC.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractPAAC)}_{a.lambda}_{a.w}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractProtFP_f = ret_extractProtFP.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractProtFP)}_{a.pc}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractProtFPGap_f = ret_extractProtFPGap.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractProtFPGap)}_{a.pc}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractQSO_f = ret_extractQSO.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractQSO)}_{a.nlag}_{a.w}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractSOCN_f = ret_extractSOCN.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractSOCN)}_{a.nlag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractScales_f = ret_extractScales.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractScales)}_{a.pc}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    var ret_extractScalesGap_f = ret_extractScalesGap.Select(a => new feature_info()
                    {
                        alphabet = alphabet_name,
                        dimension = 1,
                        category = $@"{nameof(r_protr)}",
                        source = source,
                        @group = $@"{nameof(r_protr)}_{nameof(extractScalesGap)}_{a.pc}_{a.lag}_{alphabet_name}".Replace(".", "_", StringComparison.InvariantCulture),
                        member = $@"{a.name}".Replace(".", "_", StringComparison.InvariantCulture),
                        perspective = $@"default",
                        feature_value = descriptive_stats.fix_double(a.value)
                    }).ToList();

                    // todo: replace '.'s in names

                    features.AddRange(ret_extractAPAAC_f);
                    features.AddRange(ret_extractBLOSUM_f);
                    features.AddRange(ret_extractCTDC_f);
                    features.AddRange(ret_extractCTDD_f);
                    features.AddRange(ret_extractCTDT_f);
                    features.AddRange(ret_extractCTriad_f);
                    features.AddRange(ret_extractCTriadClass_f);
                    features.AddRange(ret_extractDC_f);
                    features.AddRange(ret_extractDescScales_f);
                    features.AddRange(ret_extractFAScales_f);
                    features.AddRange(ret_extractGeary_f);
                    features.AddRange(ret_extractMDSScales_f);
                    features.AddRange(ret_extractMoran_f);
                    features.AddRange(ret_extractMoreauBroto_f);
                    features.AddRange(ret_extractPAAC_f);
                    features.AddRange(ret_extractProtFP_f);
                    features.AddRange(ret_extractProtFPGap_f);
                    features.AddRange(ret_extractQSO_f);
                    features.AddRange(ret_extractSOCN_f);
                    features.AddRange(ret_extractScales_f);
                    features.AddRange(ret_extractScalesGap_f);

                    //features.ForEach(a=>//io_proxy.WriteLine(a.ToString()));

                    if (_template_get_values == null)
                    {
                        var template = features.Select(a => new feature_info(a) { alphabet = "", source = "", feature_value = 0 }).ToList();
                        _template_get_values = template;
                    }

                    return features;
                }
            }
        }


        private static List<(string name, int lambda, double value)> template_extractAPAAC = null;

        public static List<(string name, int lambda, double value)> extractAPAAC(REngine engine, string x, int lambda_first = 1, int lambda_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_first), lambda_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_last), lambda_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractAPAAC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (template_extractAPAAC == null && (list != null && list.Count > 0))
                {
                    template_extractAPAAC = list.Select(a => (name: a.name, lambda: a.lambda, value: 0d)).ToList();
                }
                else if (template_extractAPAAC != null && (list == null || list.Count == 0))
                {
                    list = template_extractAPAAC;
                }

                return list;
            }
        }


        private static List<(string name, string submat, int k, int lag, double value)> template_extractBLOSUM = null;

        public static List<(string name, string submat, int k, int lag, double value)> extractBLOSUM(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractBLOSUM)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                var list = new List<(string name, string submat, int k, int lag, double value)>();
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

                            list.AddRange(names.Select((a, i) => (name: a, submat: submat, k: k, lag: lag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                        }
                    }
                }


                if (template_extractBLOSUM == null && (list != null && list.Count > 0))
                {
                    template_extractBLOSUM = list.Select(a => (name: a.name, submat: a.submat, k: a.k, lag: a.lag, value: 0d)).ToList();
                }
                else if (template_extractBLOSUM != null && (list == null || list.Count == 0))
                {
                    list = template_extractBLOSUM;
                }


                return list;
            }
        }


        private static List<(string name, double value)> template_extractCTDC = null;

        public static List<(string name, double value)> extractCTDC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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


                if (template_extractCTDC == null && (list != null && list.Count > 0))
                {
                    template_extractCTDC = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (template_extractCTDC != null && (list == null || list.Count == 0))
                {
                    list = template_extractCTDC;
                }

                return list;
            }
        }


        private static List<(string name, double value)> temlate_extractCTDD = null;
        public static List<(string name, double value)> extractCTDD(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDD)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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



                if (temlate_extractCTDD == null && (list != null && list.Count > 0))
                {
                    temlate_extractCTDD = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (temlate_extractCTDD != null && (list == null || list.Count == 0))
                {
                    list = temlate_extractCTDD;
                }


                return list;
            }
        }


        private static List<(string name, double value)> temlate_extractCTDT = null;

        public static List<(string name, double value)> extractCTDT(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTDT)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (temlate_extractCTDT == null && (list != null && list.Count > 0))
                {
                    temlate_extractCTDT = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (temlate_extractCTDT != null && (list == null || list.Count == 0))
                {
                    list = temlate_extractCTDT;
                }

                return list;
            }
        }


        private static List<(string name, double value)> temlate_extractCTriad = null;

        public static List<(string name, double value)> extractCTriad(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTriad)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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


                if (temlate_extractCTriad == null && (list != null && list.Count > 0))
                {
                    temlate_extractCTriad = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (temlate_extractCTriad != null && (list == null || list.Count == 0))
                {
                    list = temlate_extractCTriad;
                }

                return list;
            }
        }

        private static List<(string name, double value)> temlate_extractCTriadClass = null;


        public static List<(string name, double value)> extractCTriadClass(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractCTriadClass)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (temlate_extractCTriadClass == null && (list != null && list.Count > 0))
                {
                    temlate_extractCTriadClass = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (temlate_extractCTriadClass != null && (list == null || list.Count == 0))
                {
                    list = temlate_extractCTriadClass;
                }

                return list;
            }
        }

        private static List<(string name, double value)> temlate_extractDC = null;

        public static List<(string name, double value)> extractDC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractDC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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


                if (temlate_extractDC == null && (list != null && list.Count > 0))
                {
                    temlate_extractDC = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (temlate_extractDC != null && (list == null || list.Count == 0))
                {
                    list = temlate_extractDC;
                }

                return list;
            }
        }


        private static List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)> template_extractDescScales = null;

        public static List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)> extractDescScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractDescScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                var list = new List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)>();

                foreach (var propmat in propmats)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        for (var lag = lag_first; lag <= lag_last; lag++)
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

                            var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            for (var i = 2; i <= 4; i++)
                            {
                                var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                                var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                                var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                                pca_list.AddRange(pca_row_values.Select((a, j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                            }

                            if (names.Length != values.Length) throw new Exception();

                            list.AddRange(names.Select((a, i) => (name: a, propmat: propmat, pc: pc, lag: lag, pca: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                        }
                    }
                }


                if (template_extractDescScales == null && (list != null && list.Count > 0))
                {
                    template_extractDescScales = list.Select(a => (name: a.name, propmat: a.propmat, pc: a.pc, lag: a.lag, pca: a.pca.Select(b => (b.row_name, b.col_name, 0d)).ToList(), value: 0d)).ToList();
                }
                else if (template_extractDescScales != null && (list == null || list.Count == 0))
                {
                    list = template_extractDescScales;
                }

                return list;
            }
        }


        private static List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)> template_extractFAScales = null;

        public static List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)> extractFAScales(REngine engine, string x, int lag_first = 1, int lag_last = 2, int factors_first = 5, int factors_last = 5)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(factors_first), factors_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(factors_last), factors_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractFAScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var silent = "FALSE";
                var scale = "TRUE";

                var list = new List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)>();

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
                        var factors_rows = values_output.Skip(23).Take(3).Select(a => a.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)).ToList();
                        var factors_row_labels = factors_rows.Select(a => string.Join(" ", a.Take(a.Length - factors_col_labels.Length).ToList())).ToList();
                        var factors_row_values = factors_rows.Select(a => a.Skip(a.Length - factors_col_labels.Length).Select(b => double.Parse(b, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList()).ToList();

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
                        var p_value = double.Parse(values_output[29].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Last(), NumberStyles.Float, CultureInfo.InvariantCulture);

                        list.AddRange(names.Select((a, i) => (name: a, factors: factors, lag: lag, factors_list: factors_list, chi_sq: chi_sq, p_value: p_value, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                if (template_extractFAScales == null && (list != null && list.Count > 0))
                {
                    template_extractFAScales = list.Select(a => (name: a.name, factors: a.factors, lag: a.lag, factors_list: a.factors_list.Select(b => (row_name: b.row_name, col_name: b.col_name, value: 0d)).ToList(), chi_sq: 0d, p_value: 0d, value: 0d)).ToList();
                }
                else if (template_extractFAScales != null && (list == null || list.Count == 0))
                {
                    list = template_extractFAScales;
                }

                return list;
            }
        }



        private static List<(string name, int nlag, double value)> _template_extractGeary = null;

        public static List<(string name, int nlag, double value)> extractGeary(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractGeary)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (_template_extractGeary == null && (list != null && list.Count > 0))
                {
                    _template_extractGeary = list.Select(a => (name: a.name, nlag: a.nlag, value: 0d)).ToList();
                }
                else if (_template_extractGeary != null && (list == null || list.Count == 0))
                {
                    list = _template_extractGeary;
                }

                return list;
            }
        }


        private static List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)> _template_extractMDSScales = null;

        public static List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)> extractMDSScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractMDSScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                        var scaling_eigenvalues = values_output.Skip(1).SelectMany(a => a.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Skip(1).Select(b => double.Parse(b, NumberStyles.Float, CultureInfo.InvariantCulture)).ToArray()).ToArray();

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, k: k, lag: lag, scaling_eigenvalues: scaling_eigenvalues, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                if (_template_extractMDSScales == null && (list != null && list.Count > 0))
                {
                    _template_extractMDSScales = list.Select(a => (name: a.name, k: a.k, lag: a.lag, scaling_eigenvalues: a.scaling_eigenvalues.Select(b => 0d).ToArray(), value: 0d)).ToList();
                }
                else if (_template_extractMDSScales != null && (list == null || list.Count == 0))
                {
                    list = _template_extractMDSScales;
                }

                return list;
            }
        }


        private static List<(string name, int nlag, double value)> _template_extractMoran = null;

        public static List<(string name, int nlag, double value)> extractMoran(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractMoran)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (_template_extractMoran == null && (list != null && list.Count > 0))
                {
                    _template_extractMoran = list.Select(a => (name: a.name, nlag: a.nlag, value: 0d)).ToList();
                }
                else if (_template_extractMoran != null && (list == null || list.Count == 0))
                {
                    list = _template_extractMoran;
                }

                return list;
            }
        }



        private static List<(string name, int nlag, double value)> _template_extractMoreauBroto = null;

        public static List<(string name, int nlag, double value)> extractMoreauBroto(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractMoreauBroto)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
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

                if (_template_extractMoreauBroto == null && (list != null && list.Count > 0))
                {
                    _template_extractMoreauBroto = list.Select(a => (name: a.name, nlag: a.nlag, value: 0d)).ToList();
                }
                else if (_template_extractMoreauBroto != null && (list == null || list.Count == 0))
                {
                    list = _template_extractMoreauBroto;
                }

                return list;
            }
        }



        private static List<(string name, int lambda, double w, double value)> _template_extractPAAC = null;

        public static List<(string name, int lambda, double w, double value)> extractPAAC(REngine engine, string x, int lambda_first = 1, int lamda_last = 2, double w = 0.5)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_first), lambda_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lamda_last), lamda_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(w), w.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractPAAC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lambda, double w, double value)>();

                for (var lambda = lambda_first; lambda <= lamda_last; lambda++)
                {
                    var fn_name = nameof(extractPAAC);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = "";

                    var props = $@"c(""Hydrophobicity"", ""Hydrophilicity"", ""SideChainMass"")";
                    var customprops = "NULL";

                    var cmd = $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(lambda)} = {lambda},  {nameof(w)} = {w}, {nameof(customprops)} = {customprops})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) => (name: a, lambda: lambda, w: w, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                if (_template_extractPAAC == null && (list != null && list.Count > 0))
                {
                    _template_extractPAAC = list.Select(a => (name: a.name, lambda: a.lambda, w: a.w, value: 0d)).ToList();
                }
                else if (_template_extractPAAC != null && (list == null || list.Count == 0))
                {
                    list = _template_extractPAAC;
                }

                return list;
            }
        }



        public static void extractPSSM(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractPSSM)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            throw new NotImplementedException();
        }



        public static void extractPSSMAcc(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractPSSMAcc)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            throw new NotImplementedException();
        }



        public static void extractPSSMFeature(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractPSSMFeature)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            throw new NotImplementedException();
        }





        private static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractProtFP = null;

        public static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractProtFP(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractProtFP)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();

                var index = "c(160:165, 258:296)";
                var scale = "TRUE";
                var silent = "FALSE";

                var pc_first = 5;
                var pc_last = 5;


                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        var fn_name = nameof(extractProtFP);
                        var vr_name = $"{fn_name}_v{Key}";
                        var vr_name_output = $"{fn_name}_v{Key}";
                        var ar_ix = "";


                        var cmd = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(index)} = {index}, {nameof(pc)} = {pc},  {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";
                        var evaluate = engine.Evaluate(cmd);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacter();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_output})");

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                            pca_list.AddRange(pca_row_values.Select((a, j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                        }

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, lag: lag, pc: pc, pca_list: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }


                if (_template_extractProtFP == null && (list != null && list.Count > 0))
                {
                    _template_extractProtFP = list.Select(a => (name: a.name, lag: a.lag, pc: a.pc, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();
                }
                else if (_template_extractProtFP != null && (list == null || list.Count == 0))
                {
                    list = _template_extractProtFP;
                }

                return list;
            }
        }



        private static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractProtFPGap = null;

        public static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractProtFPGap(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractProtFPGap)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();

                var index = "c(160:165, 258:296)";
                var scale = "TRUE";
                var silent = "FALSE";

                var pc_first = 5;
                var pc_last = 5;


                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        var fn_name = nameof(extractProtFPGap);
                        var vr_name = $"{fn_name}_v{Key}";
                        var vr_name_output = $"{fn_name}_v{Key}";
                        var ar_ix = "";


                        var cmd = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(index)} = {index}, {nameof(pc)} = {pc},  {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";
                        var evaluate = engine.Evaluate(cmd);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacter();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_output})");

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                            pca_list.AddRange(pca_row_values.Select((a, j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                        }

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, lag: lag, pc: pc, pca_list: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                if (_template_extractProtFPGap == null && (list != null && list.Count > 0))
                {
                    _template_extractProtFPGap = list.Select(a => (name: a.name, lag: a.lag, pc: a.pc, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();
                }
                else if (_template_extractProtFPGap != null && (list == null || list.Count == 0))
                {
                    list = _template_extractProtFPGap;
                }

                return list;
            }
        }



        private static List<(string name, int nlag, double w, double value)> _template_extractQSO = null;

        public static List<(string name, int nlag, double w, double value)> extractQSO(REngine engine, string x, int nlag_first = 1, int nlag_last = 2, double w = 0.1)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(w), w.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractQSO)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double w, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    var fn_name = nameof(extractQSO);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = "";

                    var cmd = $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(nlag)} = {nlag},  {nameof(w)} = {w})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) => (name: a, nlag: nlag, w: w, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }


                if (_template_extractQSO == null && (list != null && list.Count > 0))
                {
                    _template_extractQSO = list.Select(a => (name: a.name, nlag: a.nlag, w: a.w, value: 0d)).ToList();
                }
                else if (_template_extractQSO != null && (list == null || list.Count == 0))
                {
                    list = _template_extractQSO;
                }

                return list;
            }
        }



        private static List<(string name, int nlag, double value)> _template_extractSOCN = null;

        public static List<(string name, int nlag, double value)> extractSOCN(REngine engine, string x, int nlag_first = 1, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractSOCN)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    var fn_name = nameof(extractSOCN);
                    var vr_name = $"{fn_name}_v{Key}";
                    var ar_ix = "";

                    var cmd = $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(nlag)} = {nlag})";
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                    var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                    var rm = engine.Evaluate($@"rm({vr_name})");

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) => (name: a, nlag: nlag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                }

                if (_template_extractSOCN == null && (list != null && list.Count > 0))
                {
                    _template_extractSOCN = list.Select(a => (name: a.name, nlag: a.nlag, value: 0d)).ToList();
                }
                else if (_template_extractSOCN != null && (list == null || list.Count == 0))
                {
                    list = _template_extractSOCN;
                }

                return list;
            }
        }





        private static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractScales = null;

        public static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractScales)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var scale = "TRUE";
                var silent = "FALSE";

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();


                for (var pc = pc_first; pc <= pc_last; pc++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        var fn_name = nameof(extractScales);
                        var vr_name = $"{fn_name}_v{Key}";
                        var ar_ix = ""; //"[[1]]";

                        var vr_name_ = $"{fn_name}_v{Key}";
                        var vr_name_output = $"{fn_name}_v{Key}";
                        var propmat = $"{fn_name}_v{Key}";

                        var cmd1 = $@"{propmat} <- t(na.omit(as.matrix(AAindex[, 7:26])))";
                        var cmd2 = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = {propmat}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacterMatrix();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_output})");
                        var rm3 = engine.Evaluate($@"rm({propmat})");

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                            pca_list.AddRange(pca_row_values.Select((a, j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                        }

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, pc: pc, lag: lag, pca_list: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                if (_template_extractScales == null && (list != null && list.Count > 0))
                {
                    _template_extractScales = list.Select(a => (name: a.name, pc: a.pc, lag: a.lag, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();
                }
                else if (_template_extractScales != null && (list == null || list.Count == 0))
                {
                    list = _template_extractScales;
                }

                return list;
            }
        }



        private static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractScalesGap = null;

        public static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractScalesGap(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractScalesGap)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var scale = "TRUE";
                var silent = "FALSE";

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();


                for (var pc = pc_first; pc <= pc_last; pc++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        var fn_name = nameof(extractScalesGap);
                        var vr_name = $"{fn_name}_v{Key}";
                        var ar_ix = ""; //"[[1]]";

                        var vr_name_ = $"{fn_name}_v{Key}";
                        var vr_name_output = $"{fn_name}_v{Key}";
                        var propmat = $"{fn_name}_v{Key}";

                        var cmd1 = $@"{propmat} <- t(na.omit(as.matrix(AAindex[, 7:26])))";
                        var cmd2 = $@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = {propmat}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )";

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                        //var dimnames = engine.Evaluate($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                        var values_output = engine.Evaluate($@"{vr_name_output}{ar_ix}").AsCharacterMatrix();
                        var rm1 = engine.Evaluate($@"rm({vr_name})");
                        var rm2 = engine.Evaluate($@"rm({vr_name_output})");
                        var rm3 = engine.Evaluate($@"rm({propmat})");

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(" ", pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();

                            pca_list.AddRange(pca_row_values.Select((a, j) => (row_name: pca_row_label, col_name: pca_col_labels[j], pca_value: a)).ToList());
                        }

                        if (names.Length != values.Length) throw new Exception();

                        list.AddRange(names.Select((a, i) => (name: a, pc: pc, lag: lag, pca_list: pca_list, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
                    }
                }

                if (_template_extractScalesGap == null && (list != null && list.Count > 0))
                {
                    _template_extractScalesGap = list.Select(a => (name: a.name, pc: a.pc, lag: a.lag, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();
                }
                else if (_template_extractScalesGap != null && (list == null || list.Count == 0))
                {
                    list = _template_extractScalesGap;
                }


                return list;
            }
        }




        private static List<(string name, double value)> template_extractTC = null;

        public static List<(string name, double value)> extractTC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? ""),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine($@"{nameof(r_protr)}.{nameof(extractTC)}({string.Join(", ", args.Select(a => $"{a.key} = \"{a.value}\"").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                var fn_name = nameof(extractTC);
                var vr_name = $"{fn_name}_v{Key}";
                var ar_ix = "";

                var cmd = $@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")";
                var evaluate = engine.Evaluate(cmd);
                var names = engine.Evaluate($@"names({vr_name}{ar_ix})").AsCharacter();
                var values = engine.Evaluate($@"{vr_name}{ar_ix}").AsNumeric();
                var rm = engine.Evaluate($@"rm({vr_name})");

                if (names.Length != values.Length) throw new Exception();

                list.AddRange(names.Select((a, i) => (name: a, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());

                if (template_extractTC == null && (list != null && list.Count > 0))
                {
                    template_extractTC = list.Select(a => (name: a.name, value: 0d)).ToList();
                }
                else if (template_extractTC != null && (list == null || list.Count == 0))
                {
                    list = template_extractTC;
                }

                return list;
            }
        }




    }
}
