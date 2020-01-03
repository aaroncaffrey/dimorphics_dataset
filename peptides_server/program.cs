﻿using System;
using dimorphics_dataset;

namespace peptides_server
{
    public static class program
    {

        public static void Test()
        {
            var x = r_peptides.get_values("AAAALLLLYYYYYY", 1);

            foreach (var a in x)
            {
                Console.WriteLine(a);
            }

            Console.ReadKey();
        }

        public static void Main(string[] args)
        {
            // tested and working with R 3.4.4

            //Test();
            //return;

            var id = int.Parse(args[0]);
            var seq = args[1];

            var x = r_peptides.get_values(seq,id);

            var y = new subsequence_classification_data.feature_info_container {feautre_info_list = x};

            var z = subsequence_classification_data.feature_info_container.serialise_json(y);

            Console.Write(z);

            //Console.ReadLine();
        }
    }
}
