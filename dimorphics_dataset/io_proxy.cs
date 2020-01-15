using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class io_proxy
    {
        public static void WriteLine(string text = "", string module = "", string func = "")
        {
            try
            {
                Console.WriteLine($@"{DateTime.Now:G} {module}.{func} -> {text}");
            }
            catch (Exception)
            {

            }
        }

        public static void CreateDirectory(string filename)//, string module_name, string function_name)
        {
            //program.WriteLine($"{module_name}.{function_name} -> {nameof(program)}.{nameof(CreateDirectory)} ( {filename} )");

            var dir = Path.GetDirectoryName(filename);

            if (!String.IsNullOrWhiteSpace(dir))
            {
                Directory.CreateDirectory(dir);
            }
            else
            {
                throw new Exception();
            }
        }

        public static string[] ReadAllLines(string filename)
        {
            while (true)
            {
                try
                {
                    var ret = File.ReadAllLines(filename);

                    return ret;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString());
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        public static void WriteAllLines(string filename, IEnumerable<string> lines, string module_name, string function_name)
        {
            //program.WriteLine($"{module_name}.{function_name} -> {nameof(program)}.{nameof(WriteAllLines)} ( {filename} )");

            CreateDirectory(filename);//, module_name, function_name);

            while (true)
            {
                try
                {
                    File.WriteAllLines(filename, lines);
                    return;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString(), nameof(io_proxy), nameof(WriteAllLines));
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        public static void AppendAllLines(string filename, IEnumerable<string> lines, string module_name, string function_name)
        {
            //program.WriteLine($"{module_name}.{function_name} -> {nameof(program)}.{nameof(AppendAllLines)} ( {filename} )");

            CreateDirectory(filename);//, module_name, function_name);

            while (true)
            {
                try
                {
                    File.AppendAllLines(filename, lines);
                    return;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString(), nameof(io_proxy), nameof(AppendAllLines));
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        public static void AppendAllText(string filename, string text, string module_name, string function_name)
        {
            //program.WriteLine($"{module_name}.{function_name} -> {nameof(program)}.{nameof(AppendAllText)} ( {filename} )");

            CreateDirectory(filename);//, module_name, function_name);

            while (true)
            {
                try
                {
                    File.AppendAllText(filename, text);
                    return;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString(), nameof(io_proxy), nameof(AppendAllText));
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        public static void WriteAllText(string filename, string text, string module_name, string function_name)
        {
            //program.WriteLine($"{module_name}.{function_name} -> {nameof(program)}.{nameof(WriteAllText)} ( {filename} )");

            CreateDirectory(filename);//, module_name, function_name);

            while (true)
            {
                try
                {
                    File.WriteAllText(filename, text);
                    return;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString(), nameof(io_proxy), nameof(WriteAllText));
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }
    }
}
