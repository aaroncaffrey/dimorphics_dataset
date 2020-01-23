using System;
using System.Globalization;
using System.IO;
using System.Linq;
using dimorphics_dataset;

namespace peptides_server
{
    public static class program
    {

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

                var x = r_peptides.get_values(0, "", "", peptide);

                foreach (var a in x)
                {
                    Console.WriteLine(a);
                }


                var y = new feature_info_container { feautre_info_list = x };

                var s = feature_info_container.serialise_json(y);

                var d = feature_info_container.deserialise(s);

            }

            Console.ReadKey();
        }

        public static void Main(string[] args)
        {
#if DEBUG
            //Test();
#endif
            if (args == null || args.Length == 0)
            {
                throw new ArgumentNullException(nameof(args));
            }

            var arg_index = 0;

            var call_count = int.Parse(args[arg_index++], NumberStyles.Integer, CultureInfo.InvariantCulture);
            var alphabet_name = args[arg_index++];
            var source_name = args[arg_index++];
            var sequence = args[arg_index++];

            var x = r_peptides.get_values(call_count, source_name, alphabet_name, sequence);

            var y = new feature_info_container {feautre_info_list = x};

            var z = feature_info_container.serialise_json(y);

            Console.Write(z);

            //Console.ReadLine();
        }
    }
}
