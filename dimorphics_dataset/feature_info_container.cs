using System;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace dimorphics_dataset
{
    public class feature_info_container
    {
        public const string module_name = nameof(feature_info_container);
        public List<feature_info> feautre_info_list;

        public static string serialise_json(feature_info_container feature_info_container)
        {
            var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All, };
            var serialised_serialise_json = JsonConvert.SerializeObject(feature_info_container, json_settings);
            return serialised_serialise_json;
        }

        public static feature_info_container deserialise(string serialized_json)
        {
            
            const string method_name = nameof(deserialise);

            var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All };
            feature_info_container feature_info_container;

            try
            {
                feature_info_container = JsonConvert.DeserializeObject<feature_info_container>(serialized_json, json_settings);
            }
            catch (Exception e)
            {
                feature_info_container = null;

                io_proxy.log_exception(e, /*program.string_debug*/($@""), module_name, method_name);
            }

            return feature_info_container;
        }
    }
}
