using System;
using System.Globalization;
using dimorphics_dataset;

namespace protr_server
{
    public static class program
    {
        public static void Test()
        {
            var x = r_protr.get_values("AAAALLLLYYYYYY", 1);

            foreach (var a in x)
            {
                Console.WriteLine(a);
            }


            var y = new subsequence_classification_data.feature_info_container { feautre_info_list = x };

            var s = subsequence_classification_data.feature_info_container.serialise_json(y);

            var d = subsequence_classification_data.feature_info_container.deserialise(s);

            Console.ReadKey();
        }

        public static void Main(string[] args)
        {
#if DEBUG
            args = new string[] {"0", "AAA"};
#endif
            
            var call_count = int.Parse(args[0], NumberStyles.Integer, CultureInfo.InvariantCulture);
            
            var sequence = args[1];

            var x = r_protr.get_values(sequence, call_count);

            var y = new subsequence_classification_data.feature_info_container {feautre_info_list = x};

            var z = subsequence_classification_data.feature_info_container.serialise_json(y);

            Console.Write(z);
        }
    }
}
