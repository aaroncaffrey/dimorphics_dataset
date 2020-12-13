using System;
using System.Globalization;
using dimorphics_dataset;

namespace protr_server
{
    internal static class program
    {/*
#if DEBUG
        private static Random random = new Random();
        internal static string random_peptide(int length)
        {
            const string chars = /*program.string_debug* /($@"ACDEFGHIKLMNPQRSTVWY";
            return new string(Enumerable.Repeat(chars, length).Select(s => s[random.Next(s.Length)]).ToArray());
        }

        internal static void Test()
        {
            for (var length = 3; length <= 10; length++)
            {
                var peptide = random_peptide(length);

                var x = r_protr.get_values(0, /*program.string_debug* /($@"datasource", /*program.string_debug* /($@"Overall", peptide);

                io_proxy.WriteLine(/*program.string_debug* /($@"input size: " + length + /*program.string_debug* /($@", output size: " + x.Count);
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

        internal static string string_debug(string str)
        {
#if DEBUG
            var i_str = string.IsInterned(str);
            if (i_str == null)
            {
                i_str = string.Intern(str);
                Console.WriteLine($"{nameof(string_debug)}: {i_str}");
                return i_str;
            }
#endif
            return str;
        }

        internal static void Main(string[] args)
        {
            // protr_server.program.Main

            if (args == null || args.Length == 0)
            {
                throw new ArgumentNullException(nameof(args));
            }
            
            var arg_index = 0;

            var call_count = int.Parse(args[arg_index++], NumberStyles.Integer, NumberFormatInfo.InvariantInfo);
            var alphabet_name = /*program.string_debug*/($@"");//args[arg_index++];
            var source_name = /*program.string_debug*/($@"");//args[arg_index++];
            var sequence = /*program.string_debug*/(args[arg_index++]);

            var x = r_protr.get_values(call_count, source_name, alphabet_name, sequence);

            var y = new feature_info_container {feautre_info_list = x};

            var z = feature_info_container.serialise_json(y);

            Console.Write(z);
        }
    }
}
