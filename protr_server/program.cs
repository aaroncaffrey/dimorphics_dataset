using System;
using System.Globalization;
using System.Linq;
using dimorphics_dataset;

namespace protr_server
{
    public static class program
    {/*
#if DEBUG
        private static Random random = new Random();
        public static string random_peptide(int length)
        {
            const string chars = "ACDEFGHIKLMNPQRSTVWY";
            return new string(Enumerable.Repeat(chars, length).Select(s => s[random.Next(s.Length)]).ToArray());
        }

        public static void Test()
        {
            for (var length = 3; length <= 10; length++)
            {
                var peptide = random_peptide(length);

                var x = r_protr.get_values(0, "datasource", "Overall", peptide);

                io_proxy.WriteLine("input size: " + length + ", output size: " + x.Count);
                foreach (var a in x)
                {
                    //io_proxy.WriteLine(a.ToString());
                }


                var y = new feature_info_container {feautre_info_list = x};

                var s = feature_info_container.serialise_json(y);

                var d = feature_info_container.deserialise(s);

            }

            Console.ReadKey();
        }
#endif*/

        public static void Main(string[] args)
        {
            if (args == null || args.Length == 0)
            {
                throw new ArgumentNullException(nameof(args));
            }
            
            var arg_index = 0;

            var call_count = int.Parse(args[arg_index++], NumberStyles.Integer, CultureInfo.InvariantCulture);
            var alphabet_name = "";//args[arg_index++];
            var source_name = "";//args[arg_index++];
            var sequence = args[arg_index++];

            var x = r_protr.get_values(call_count, source_name, alphabet_name, sequence);

            var y = new feature_info_container {feautre_info_list = x};

            var z = feature_info_container.serialise_json(y);

            Console.Write(z);
        }
    }
}
