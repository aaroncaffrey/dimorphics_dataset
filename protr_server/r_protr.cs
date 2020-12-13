using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using dimorphics_dataset;
using RDotNet;

namespace protr_server
{
    internal static class r_protr
    {
        internal static readonly object engine_lock = new object();
        private static uint _key = 1;
        private static readonly object _key_lock = new object();
        private static int _id;
        private static List<feature_info> _template_get_values;
        private static readonly object _get_values_lock = new object();

        internal static string Key
        {
            get
            {
                lock (_key_lock)
                {
                    _key++;
                    return /*program.string_debug*/($@"{nameof(r_protr)}_{_id}_{_key}");
                }
            }
        }

        private static bool need_init = true;

        internal static REngine init_r()
        {
            var rinit = new StartupParameter {Quiet = true, Interactive = false, RHome = /*program.string_debug*/($@"C:\Program Files\R\R-3.6.2\")};
            var engine1 = REngine.GetInstance(Path.Combine(rinit.RHome, /*program.string_debug*/($@"bin\x64\R.dll")), true, rinit);

            if (need_init)
            {
                need_init = false;

                var r_init_cmds = /*program.string_debug*/($@"
                    #install.packages(""devtools"")
                    #library(devtools)
                    #install_github(""https://github.com/nanxstats/protr"")
                    library(protr)
                ");

                r_init_cmds.Split(new char[] { '\r', '\n' }).Where(a => !string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith(/*program.string_debug*/($@"#"), StringComparison.Ordinal)).ToList().ForEach(a => engine1.Evaluate(a));
            }

            return engine1;
        }

        internal static List<feature_info> get_values(int id, string source, string alphabet_name, string sequence, int min_sequence_length = 1)
        {
            lock (_get_values_lock)
            {
                r_protr._id = id;

                var engine = init_r();
                {
                    if (string.IsNullOrWhiteSpace(sequence) || sequence.Length < min_sequence_length)
                    {
                        if (_template_get_values == null)
                        {
                            _template_get_values = get_values(id, source, alphabet_name, /*program.string_debug*/($@"ALG"));
                            _template_get_values = _template_get_values.Select(a => new feature_info(a) { alphabet = alphabet_name, source = source, feature_value = 0 }).ToList();
                        }

                        return _template_get_values;
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

                    var ret_extractAPAAC_f = ret_extractAPAAC.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractAPAAC)}_{ds_stat.lambda}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractBLOSUM_f = ret_extractBLOSUM.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractBLOSUM)}_{ds_stat.submat}_{ds_stat.k}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractCTDC_f = ret_extractCTDC.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractCTDC)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractCTDD_f = ret_extractCTDD.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractCTDD)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractCTDT_f = ret_extractCTDT.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractCTDT)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractCTriad_f = ret_extractCTriad.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractCTriad)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractCTriadClass_f = ret_extractCTriadClass.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractCTriadClass)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractDC_f = ret_extractDC.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractDC)}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractDescScales_f = ret_extractDescScales.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractDescScales)}_{ds_stat.propmat}_{ds_stat.pc}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractFAScales_f = ret_extractFAScales.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractFAScales)}_{ds_stat.factors}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractGeary_f = ret_extractGeary.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractGeary)}_{ds_stat.nlag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractMDSScales_f = ret_extractMDSScales.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractMDSScales)}_{ds_stat.k}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractMoran_f = ret_extractMoran.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractMoran)}_{ds_stat.nlag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractMoreauBroto_f = ret_extractMoreauBroto.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractMoreauBroto)}_{ds_stat.nlag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractPAAC_f = ret_extractPAAC.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractPAAC)}_{ds_stat.lambda}_{ds_stat.w}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractProtFP_f = ret_extractProtFP.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractProtFP)}_{ds_stat.pc}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractProtFPGap_f = ret_extractProtFPGap.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractProtFPGap)}_{ds_stat.pc}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractQSO_f = ret_extractQSO.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractQSO)}_{ds_stat.nlag}_{ds_stat.w}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractSOCN_f = ret_extractSOCN.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractSOCN)}_{ds_stat.nlag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractScales_f = ret_extractScales.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractScales)}_{ds_stat.pc}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
                    }).ToList();

                    var ret_extractScalesGap_f = ret_extractScalesGap.Select(ds_stat => new feature_info()
                    {
                        alphabet = alphabet_name,
                        stats = "", dimension = 1,
                        category = /*program.string_debug*/($@"{nameof(r_protr)}"),
                        source = source,
                        @group = /*program.string_debug*/($@"{nameof(r_protr)}_{nameof(extractScalesGap)}_{ds_stat.pc}_{ds_stat.lag}_{alphabet_name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        member = /*program.string_debug*/($@"{ds_stat.name}".Replace(/*program.string_debug*/($@"."), /*program.string_debug*/($@"_"), StringComparison.Ordinal)),
                        perspective = /*program.string_debug*/($@"default"),
                        feature_value = descriptive_stats.fix_double(ds_stat.value)
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
                        var template = features.Select(a => new feature_info(a) { alphabet = /*program.string_debug*/($@""), source = /*program.string_debug*/($@""), feature_value = 0 }).ToList();
                        _template_get_values = template;
                    }

                    return features;
                }
            }
        }


        private static List<(string name, int lambda, double value)> template_extractAPAAC;

        internal static  List<(string name, int lambda, double value)> extractAPAAC(REngine engine, string x, int lambda_first = 1, int lambda_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lambda_first || x.Length <= lambda_last)
            {
                if (template_extractAPAAC == null)
                {
                    var seq_len = (lambda_last > lambda_first ? lambda_last : lambda_first) + 1;
                    var temp_str = /*program.string_debug*/(string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len));
                    template_extractAPAAC = extractAPAAC(engine, temp_str, lambda_first, lambda_last);
                    template_extractAPAAC = template_extractAPAAC.Select(a => (name: a.name, lambda: a.lambda, value: 0d)).ToList();
                }

                return template_extractAPAAC;
            }
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_first), lambda_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_last), lambda_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractAPAAC)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lambda, double value)>();

                for (var lambda = lambda_first; lambda <= lambda_last; lambda++)
                {
                    var f = nameof(extractAPAAC);
                    var k = Key;
                    var v = /*program.string_debug*/($@"{f}_v{k}");
                    var ai_ix = /*program.string_debug*/($@""); //"[[1]]";

                    var evaluate = engine.Evaluate(/*program.string_debug*/($@"{v} <- {f}({nameof(x)} = ""{x}"", {nameof(lambda)} = {lambda})"));
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({v}{ai_ix})")).AsCharacter();
                    //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({v}{ai})");
                    //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({v}{ai})");
                    //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({v}{ai})");
                    var values = engine.Evaluate(/*program.string_debug*/($@"{v}{ai_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({v})"));

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


        private static List<(string name, string submat, int k, int lag, double value)> template_extractBLOSUM;

        internal static List<(string name, string submat, int k, int lag, double value)> extractBLOSUM(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (template_extractBLOSUM == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    template_extractBLOSUM = extractBLOSUM(engine, temp_str, lag_first, lag_last);
                    template_extractBLOSUM = template_extractBLOSUM.Select(a => (name: a.name, submat: a.submat, k: a.k, lag: a.lag, value: 0d)).ToList();

                }

                return template_extractBLOSUM;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractBLOSUM)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var submats = new string[]
                {
                    /*program.string_debug*/($@"AABLOSUM45"), /*program.string_debug*/($@"AABLOSUM50"), /*program.string_debug*/($@"AABLOSUM62"), /*program.string_debug*/($@"AABLOSUM80"), /*program.string_debug*/($@"AABLOSUM100"), /*program.string_debug*/($@"AAPAM30"), /*program.string_debug*/($@"AAPAM40"),
                    /*program.string_debug*/($@"AAPAM70"), /*program.string_debug*/($@"AAPAM120"), /*program.string_debug*/($@"AAPAM250")
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
                            const string fn_name = nameof(extractBLOSUM);
                            var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                            var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                            var evaluate =
                                engine.Evaluate(
                                    /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(submat)} = ""{submat}"", {nameof(k)} = {k}, {nameof(lag)} = {lag})"));
                            var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                            //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                            //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                            //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                            var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                            var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, double value)> template_extractCTDC;

        internal static List<(string name, double value)> extractCTDC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractCTDC)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractCTDC);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                var evaluate = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, double value)> temlate_extractCTDD;
        internal static List<(string name, double value)> extractCTDD(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractCTDD)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractCTDD);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                var evaluate = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, double value)> temlate_extractCTDT;

        internal static List<(string name, double value)> extractCTDT(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractCTDT)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractCTDT);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                var evaluate = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, double value)> temlate_extractCTriad;

        internal static List<(string name, double value)> extractCTriad(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractCTriad)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractCTriad);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@"");

                var evaluate = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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

        private static List<(string name, double value)> temlate_extractCTriadClass;


        internal static List<(string name, double value)> extractCTriadClass(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractCTriadClass)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractCTriadClass);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@"");

                var aaclass = /*program.string_debug*/($@"{fn_name}_v{Key}");

                var evaluate1 = engine.Evaluate(/*program.string_debug*/($@"{aaclass} <- list( c(""G"", ""A"", ""S"", ""T"", ""P"", ""D"", ""C""), c(""N"", ""V"", ""E"", ""Q"", ""I"", ""L""), c(""M"", ""H"", ""K"", ""F"", ""R"", ""Y"", ""W"") )"));
                var evaluate2 = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(aaclass)} = {aaclass})"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({aaclass})"));
                var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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

        private static List<(string name, double value)> temlate_extractDC;

        internal static List<(string name, double value)> extractDC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractDC)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractDC);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                var evaluate = engine.Evaluate(/*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")"));
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)> template_extractDescScales;

        internal static List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)> extractDescScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (template_extractDescScales == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    template_extractDescScales = extractDescScales(engine, temp_str, lag_first, lag_last);
                    template_extractDescScales = template_extractDescScales.Select(a => (name: a.name, propmat: a.propmat, pc: a.pc, lag: a.lag, pca: a.pca.Select(b => (b.row_name, b.col_name, 0d)).ToList(), value: 0d)).ToList();

                }

                return template_extractDescScales;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractDescScales)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var propmats = new string[]
                {
                    /*program.string_debug*/($@"AAMOE2D"), /*program.string_debug*/($@"AAMOE3D"), /*program.string_debug*/($@"AACPSA"), /*program.string_debug*/($@"AADescAll"), /*program.string_debug*/($@"AA2DACOR"), /*program.string_debug*/($@"AA3DMoRSE"), /*program.string_debug*/($@"AAACF"), /*program.string_debug*/($@"AABurden"), /*program.string_debug*/($@"AAConn"),
                    /*program.string_debug*/($@"AAConst"), /*program.string_debug*/($@"AAEdgeAdj"), /*program.string_debug*/($@"AAEigIdx"), /*program.string_debug*/($@"AAFGC"), /*program.string_debug*/($@"AAGeom"), /*program.string_debug*/($@"AAGETAWAY"), /*program.string_debug*/($@"AAInfo"), /*program.string_debug*/($@"AAMolProp"),
                    /*program.string_debug*/($@"AARandic"), /*program.string_debug*/($@"AARDF"), /*program.string_debug*/($@"AATopo"), /*program.string_debug*/($@"AATopoChg"), /*program.string_debug*/($@"AAWalk"), /*program.string_debug*/($@"AAWHIM")
                };

                var index = /*program.string_debug*/($@"NULL");
                var silent = /*program.string_debug*/($@"FALSE");
                var scale = /*program.string_debug*/($@"TRUE");

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string name, string propmat, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca, double value)>();

                foreach (var propmat in propmats)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        for (var lag = lag_first; lag <= lag_last; lag++)
                        {


                            const string fn_name = nameof(extractDescScales);
                            var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                            var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                            var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");

                            var cmd = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = ""{propmat}"", {nameof(index)} = {index}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");

                            var evaluate = engine.Evaluate(cmd);
                            var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                            //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                            //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                            //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                            var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                            var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacterMatrix();
                            var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                            var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));

                            var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                            var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                            for (var i = 2; i <= 4; i++)
                            {
                                var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                                var pca_row_label = string.Join(/*program.string_debug*/($@" "), pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                                var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList();

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


        private static List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)> template_extractFAScales;

        internal static List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)> extractFAScales(REngine engine, string x, int lag_first = 1, int lag_last = 2, int factors_first = 5, int factors_last = 5)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (template_extractFAScales == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    template_extractFAScales = extractFAScales(engine, temp_str, lag_first, lag_last, factors_first, factors_last);
                    template_extractFAScales = template_extractFAScales.Select(a => (name: a.name, factors: a.factors, lag: a.lag, factors_list: a.factors_list.Select(b => (row_name: b.row_name, col_name: b.col_name, value: 0d)).ToList(), chi_sq: 0d, p_value: 0d, value: 0d)).ToList();

                }

                return template_extractFAScales;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(factors_first), factors_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(factors_last), factors_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractFAScales)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var silent = /*program.string_debug*/($@"FALSE");
                var scale = /*program.string_debug*/($@"TRUE");

                var list = new List<(string name, int factors, int lag, List<(string row_name, string col_name, double value)> factors_list, double chi_sq, double p_value, double value)>();

                for (var factors = factors_first; factors <= factors_last; factors++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        const string fn_name = nameof(extractFAScales);
                        var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var ar_ix = /*program.string_debug*/($@"");

                        var vr_name_tprops = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");

                        var cmd1 = /*program.string_debug*/($@"{vr_name_tprops} <- AATopo[, c(37:41, 43:47)]"); // select a set of topological descriptors
                        var cmd2 = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}(x = ""{x}"", propmat = {vr_name_tprops}, {nameof(factors)} = {factors}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacter();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_tprops})"));
                        var rm3 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));

                        if (names.Length != values.Length) throw new Exception();

                        var uniqueness_labels = values_output[6].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var uniqueness_values = values_output[7].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                        // note: loadings_rows cannot be parsed correctly to be useful
                        var loadings_col_labels = values_output[10].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var loadings_rows = values_output.Skip(11).Take(uniqueness_labels.Length).Select(a => a.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)).ToList();

                        var factors_col_labels = values_output[22].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        var factors_rows = values_output.Skip(23).Take(3).Select(a => a.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries)).ToList();
                        var factors_row_labels = factors_rows.Select(a => string.Join(/*program.string_debug*/($@" "), a.Take(a.Length - factors_col_labels.Length).ToList())).ToList();
                        var factors_row_values = factors_rows.Select(a => a.Skip(a.Length - factors_col_labels.Length).Select(b => double.Parse(b, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList()).ToList();

                        var factors_list = new List<(string row_name, string col_name, double value)>();
                        for (var r = 0; r < factors_row_labels.Count; r++)
                        {
                            for (var c = 0; c < factors_col_labels.Length; c++)
                            {
                                factors_list.Add((factors_row_labels[r], factors_col_labels[c], factors_row_values[r][c]));
                            }
                        }
                        //var values_output2 = values_output[0,]

                        var chi_sq = double.Parse(values_output[28].Split()[5], NumberStyles.Float, NumberFormatInfo.InvariantInfo);
                        var p_value = double.Parse(values_output[29].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Last(), NumberStyles.Float, NumberFormatInfo.InvariantInfo);

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



        private static List<(string name, int nlag, double value)> _template_extractGeary;

        internal static List<(string name, int nlag, double value)> extractGeary(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractGeary)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    const string fn_name = nameof(extractGeary);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@"");

                    var props = /*program.string_debug*/($@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"", ""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")");
                    var customprops = /*program.string_debug*/($@"NULL");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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


        private static List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)> _template_extractMDSScales;

        internal static List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)> extractMDSScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
            
            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (_template_extractMDSScales == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractMDSScales = extractMDSScales(engine, temp_str, lag_first, lag_last);
                    _template_extractMDSScales = _template_extractMDSScales.Select(a => (name: a.name, k: a.k, lag: a.lag, scaling_eigenvalues: a.scaling_eigenvalues.Select(b => 0d).ToArray(), value: 0d)).ToList();
                }

                return _template_extractMDSScales;
            }
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractMDSScales)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int k, int lag, double[] scaling_eigenvalues, double value)>();

                const string fn_name = nameof(extractMDSScales);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@"");

                var vr_name_tprops = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");

                var scale = /*program.string_debug*/($@"TRUE");
                var silent = /*program.string_debug*/($@"FALSE");

                var k_first = 5;
                var k_last = 5;

                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var k = k_first; k <= k_last; k++)
                    {
                        var cmd1 = /*program.string_debug*/($@"{vr_name_tprops} <- AATopo[, c(37:41, 43:47)]");
                        var cmd2 = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", propmat = {vr_name_tprops}, {nameof(k)} = {k}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");
                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacter();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_tprops})"));
                        var rm3 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));

                        var scaling_eigenvalues = values_output.Skip(1).SelectMany(a => a.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Skip(1).Select(b => double.Parse(b, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToArray()).ToArray();

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


        private static List<(string name, int nlag, double value)> _template_extractMoran;

        internal static List<(string name, int nlag, double value)> extractMoran(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractMoran)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    const string fn_name = nameof(extractMoran);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                    var props = /*program.string_debug*/($@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"", ""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")");
                    var customprops = /*program.string_debug*/($@"NULL");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

                    if (names.Length != values.Length) throw new Exception();

                    list.AddRange(names.Select((a, i) => (name: a, nlag: nlag, value: double.IsNaN(values[i]) ? 0 : values[i])).ToList());
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



        private static List<(string name, int nlag, double value)> _template_extractMoreauBroto;

        internal static List<(string name, int nlag, double value)> extractMoreauBroto(REngine engine, string x, int nlag_first = 2, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractMoreauBroto)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    const string fn_name = nameof(extractMoran);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@""); //"[[1]]";


                    var props = /*program.string_debug*/($@"c(""CIDH920105"", ""BHAR880101"", ""CHAM820101"", ""CHAM820102"", ""CHOC760101"", ""BIGC670101"", ""CHAM810101"", ""DAYM780201"")");
                    var customprops = /*program.string_debug*/($@"NULL");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(nlag)} = {nlag}, {nameof(customprops)} = {customprops})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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



        private static List<(string name, int lambda, double w, double value)> _template_extractPAAC;

        internal static List<(string name, int lambda, double w, double value)> extractPAAC(REngine engine, string x, int lambda_first = 1, int lamda_last = 2, double w = 0.5)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lambda_first || x.Length <= lamda_last)
            {
                if (_template_extractPAAC == null)
                {
                    var seq_len = (lamda_last > lambda_first ? lamda_last : lambda_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractPAAC = extractPAAC(engine, temp_str, lambda_first, lamda_last);
                    _template_extractPAAC = _template_extractPAAC.Select(a => (name: a.name, lambda: a.lambda, w: a.w, value: 0d)).ToList();
                }

                return _template_extractPAAC;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lambda_first), lambda_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lamda_last), lamda_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(w), w.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractPAAC)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lambda, double w, double value)>();

                for (var lambda = lambda_first; lambda <= lamda_last; lambda++)
                {
                    const string fn_name = nameof(extractPAAC);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@"");

                    var props = /*program.string_debug*/($@"c(""Hydrophobicity"", ""Hydrophilicity"", ""SideChainMass"")");
                    var customprops = /*program.string_debug*/($@"NULL");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(props)} = {props}, {nameof(lambda)} = {lambda},  {nameof(w)} = {w}, {nameof(customprops)} = {customprops})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                    //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                    //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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



        internal static void extractPSSM(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractPSSM)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            throw new NotImplementedException();
        }



        internal static void extractPSSMAcc(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractPSSMAcc)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            throw new NotImplementedException();
        }



        internal static void extractPSSMFeature(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return;

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractPSSMFeature)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            throw new NotImplementedException();
        }





        private static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractProtFP;

        internal static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractProtFP(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (_template_extractProtFP == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractProtFP = extractProtFP(engine, temp_str, lag_first, lag_last);
                    _template_extractProtFP = _template_extractProtFP.Select(a => (name: a.name, lag: a.lag, pc: a.pc, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();
                }

                return _template_extractProtFP;
            }
                 
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractProtFP)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();

                var index = /*program.string_debug*/($@"c(160:165, 258:296)");
                var scale = /*program.string_debug*/($@"TRUE");
                var silent = /*program.string_debug*/($@"FALSE");

                var pc_first = 5;
                var pc_last = 5;


                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        const string fn_name = nameof(extractProtFP);
                        var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var ar_ix = /*program.string_debug*/($@"");


                        var cmd = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(index)} = {index}, {nameof(pc)} = {pc},  {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");
                        var evaluate = engine.Evaluate(cmd);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacter();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(/*program.string_debug*/($@" "), pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList();

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



        private static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractProtFPGap;

        internal static List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractProtFPGap(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (_template_extractProtFPGap == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractProtFPGap = extractProtFPGap(engine, temp_str, lag_first, lag_last);
                    _template_extractProtFPGap = _template_extractProtFPGap.Select(a => (name: a.name, lag: a.lag, pc: a.pc, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();

                }

                return _template_extractProtFPGap;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractProtFPGap)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int lag, int pc, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();

                var index = /*program.string_debug*/($@"c(160:165, 258:296)");
                var scale = /*program.string_debug*/($@"TRUE");
                var silent = /*program.string_debug*/($@"FALSE");

                var pc_first = 5;
                var pc_last = 5;


                for (var lag = lag_first; lag <= lag_last; lag++)
                {
                    for (var pc = pc_first; pc <= pc_last; pc++)
                    {
                        const string fn_name = nameof(extractProtFPGap);
                        var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var ar_ix = /*program.string_debug*/($@"");


                        var cmd = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(index)} = {index}, {nameof(pc)} = {pc},  {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");
                        var evaluate = engine.Evaluate(cmd);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacter();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(/*program.string_debug*/($@" "), pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList();

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



        private static List<(string name, int nlag, double w, double value)> _template_extractQSO;

        internal static List<(string name, int nlag, double w, double value)> extractQSO(REngine engine, string x, int nlag_first = 1, int nlag_last = 2, double w = 0.1)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= nlag_first || x.Length <= nlag_last)
            {
                if (_template_extractQSO == null)
                {
                    var seq_len = (nlag_last > nlag_first ? nlag_last : nlag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractQSO = extractQSO(engine, temp_str, nlag_first, nlag_last, w);
                    _template_extractQSO = _template_extractQSO.Select(a => (name: a.name, nlag: a.nlag, w: a.w, value: 0d)).ToList();
                }

                return _template_extractQSO;
            }
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
                (nameof(w), w.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractQSO)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double w, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    const string fn_name = nameof(extractQSO);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@"");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(nlag)} = {nlag},  {nameof(w)} = {w})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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



        private static List<(string name, int nlag, double value)> _template_extractSOCN;

        internal static List<(string name, int nlag, double value)> extractSOCN(REngine engine, string x, int nlag_first = 1, int nlag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= nlag_first || x.Length <= nlag_last)
            {
                if (_template_extractSOCN == null)
                {
                    var seq_len = (nlag_last > nlag_first ? nlag_last : nlag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractSOCN = extractSOCN(engine, temp_str, nlag_first, nlag_last);
                    _template_extractSOCN = _template_extractSOCN.Select(a => (name: a.name, nlag: a.nlag, value: 0d)).ToList();
                }

                return _template_extractSOCN;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_first), nlag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(nlag_last), nlag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractSOCN)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, int nlag, double value)>();

                for (var nlag = nlag_first; nlag <= nlag_last; nlag++)
                {
                    const string fn_name = nameof(extractSOCN);
                    var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                    var ar_ix = /*program.string_debug*/($@"");

                    var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(nlag)} = {nlag})");
                    var evaluate = engine.Evaluate(cmd);
                    var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                    var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                    var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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





        private static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractScales;

        internal static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractScales(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;

            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (_template_extractScales == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractScales = extractScales(engine, temp_str, lag_first, lag_last);
                    _template_extractScales = _template_extractScales.Select(a => (name: a.name, pc: a.pc, lag: a.lag, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();

                }

                return _template_extractScales;
            }
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractScales)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var scale = /*program.string_debug*/($@"TRUE");
                var silent = /*program.string_debug*/($@"FALSE");

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();


                for (var pc = pc_first; pc <= pc_last; pc++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        const string fn_name = nameof(extractScales);
                        var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                        var vr_name_ = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var propmat = /*program.string_debug*/($@"{fn_name}_v{Key}");

                        var cmd1 = /*program.string_debug*/($@"{propmat} <- t(na.omit(as.matrix(AAindex[, 7:26])))");
                        var cmd2 = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = {propmat}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacterMatrix();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));
                        var rm3 = engine.Evaluate(/*program.string_debug*/($@"rm({propmat})"));

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(/*program.string_debug*/($@" "), pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList();

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



        private static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> _template_extractScalesGap;

        internal static List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)> extractScalesGap(REngine engine, string x, int lag_first = 1, int lag_last = 2)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
            
            if (x.Length <= lag_first || x.Length <= lag_last)
            {
                if (_template_extractScalesGap == null)
                {
                    var seq_len = (lag_last > lag_first ? lag_last : lag_first) + 1;
                    var temp_str = string.Join(/*program.string_debug*/($@""), Enumerable.Repeat(/*program.string_debug*/($@"ALG"), seq_len).ToList()).Substring(0, seq_len);
                    _template_extractScalesGap = extractScalesGap(engine, temp_str, lag_first, lag_last);
                    _template_extractScalesGap = _template_extractScalesGap.Select(a => (name: a.name, pc: a.pc, lag: a.lag, pca_list: a.pca_list.Select(b => (row_name: b.row_name, col_name: b.col_name, pca_value: 0d)).ToList(), value: 0d)).ToList();

                }

                return _template_extractScalesGap;
            }

#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_first), lag_first.ToString(CultureInfo.InvariantCulture)),
                (nameof(lag_last), lag_last.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractScalesGap)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var scale = /*program.string_debug*/($@"TRUE");
                var silent = /*program.string_debug*/($@"FALSE");

                var pc_first = 5;
                var pc_last = 5;

                var list = new List<(string name, int pc, int lag, List<(string row_name, string col_name, double pca_value)> pca_list, double value)>();


                for (var pc = pc_first; pc <= pc_last; pc++)
                {
                    for (var lag = lag_first; lag <= lag_last; lag++)
                    {
                        const string fn_name = nameof(extractScalesGap);
                        var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var ar_ix = /*program.string_debug*/($@""); //"[[1]]";

                        var vr_name_ = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var vr_name_output = /*program.string_debug*/($@"{fn_name}_v{Key}");
                        var propmat = /*program.string_debug*/($@"{fn_name}_v{Key}");

                        var cmd1 = /*program.string_debug*/($@"{propmat} <- t(na.omit(as.matrix(AAindex[, 7:26])))");
                        var cmd2 = /*program.string_debug*/($@"{vr_name_output} <- capture.output( {vr_name} <- {fn_name}({nameof(x)} = ""{x}"", {nameof(propmat)} = {propmat}, {nameof(pc)} = {pc}, {nameof(lag)} = {lag}, {nameof(scale)} = {scale}, {nameof(silent)} = {silent}) )");

                        var evaluate1 = engine.Evaluate(cmd1);
                        var evaluate2 = engine.Evaluate(cmd2);
                        var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                        //var dimnames = engine.Evaluate(/*program.string_debug*/($@"dimnames({vr_name}{ar_ix})");
                        //var rownames = engine.Evaluate(/*program.string_debug*/($@"rownames({vr_name}{ar_ix})");
                        //var colnames = engine.Evaluate(/*program.string_debug*/($@"colnames({vr_name}{ar_ix})");
                        var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                        var values_output = engine.Evaluate(/*program.string_debug*/($@"{vr_name_output}{ar_ix}")).AsCharacterMatrix();
                        var rm1 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));
                        var rm2 = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name_output})"));
                        var rm3 = engine.Evaluate(/*program.string_debug*/($@"rm({propmat})"));

                        var pca_list = new List<(string row_name, string col_name, double pca_value)>();

                        var pca_col_labels = values_output[1, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        for (var i = 2; i <= 4; i++)
                        {
                            var pca_val = values_output[i, 0].Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                            var pca_row_label = string.Join(/*program.string_debug*/($@" "), pca_val.Take(pca_val.Length - pca_col_labels.Length).ToArray());
                            var pca_row_values = pca_val.Skip(pca_val.Length - pca_col_labels.Length).Select(a => double.Parse(a, NumberStyles.Float, NumberFormatInfo.InvariantInfo)).ToList();

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




        private static List<(string name, double value)> template_extractTC;

        internal static List<(string name, double value)> extractTC(REngine engine, string x)
        {
            if (engine == null || string.IsNullOrWhiteSpace(x)) return default;
#if DEBUG
            var args = new List<(string key, string value)>()
            {
                (nameof(engine), engine?.ToString() ?? /*program.string_debug*/($@"")),
                (nameof(x), x.ToString(CultureInfo.InvariantCulture)),
            };
            //io_proxy.WriteLine(/*program.string_debug*/($@"{nameof(r_protr)}.{nameof(extractTC)}({string.Join(/*program.string_debug*/($@", "), args.Select(a => /*program.string_debug*/($@"{a.key} = ""{a.value}""").ToList())})");
#endif
            lock (engine_lock)
            {
                var list = new List<(string name, double value)>();

                const string fn_name = nameof(extractTC);
                var vr_name = /*program.string_debug*/($@"{fn_name}_v{Key}");
                var ar_ix = /*program.string_debug*/($@"");

                var cmd = /*program.string_debug*/($@"{vr_name} <- {fn_name}({nameof(x)} = ""{x}"")");
                var evaluate = engine.Evaluate(cmd);
                var names = engine.Evaluate(/*program.string_debug*/($@"names({vr_name}{ar_ix})")).AsCharacter();
                var values = engine.Evaluate(/*program.string_debug*/($@"{vr_name}{ar_ix}")).AsNumeric();
                var rm = engine.Evaluate(/*program.string_debug*/($@"rm({vr_name})"));

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
