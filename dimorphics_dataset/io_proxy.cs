using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    public static class io_proxy
    {
        private static readonly object _console_lock = new object();

        public static void WriteLine(string text = "", string module_name = "", string function_name = "", bool use_lock = false)
        {
            if (!program.verbose) return;

            try
            {
                var pid = Process.GetCurrentProcess().Id;
                var thread_id = Thread.CurrentThread.ManagedThreadId;
                var task_id = Task.CurrentId ?? 0;

                if (use_lock)
                {
                    lock (_console_lock)
                    {
                        Console.WriteLine($@"{DateTime.Now:G} {pid:000000}.{thread_id:000000}.{task_id:000000} {module_name}.{function_name} -> {text}");
                    }
                }
                else
                {
                    Console.WriteLine($@"{DateTime.Now:G} {pid:000000}.{thread_id:000000}.{task_id:000000} {module_name}.{function_name} -> {text}");
                }
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

        public static string[] ReadAllLines(string filename, string module_name="", string function_name="")
        {
            io_proxy.WriteLine($"{module_name}.{function_name} -> ( {filename} )", nameof(io_proxy), nameof(ReadAllLines));

            while (true)
            {
                try
                {
                    var ret = File.ReadAllLines(filename);

                    return ret;
                }
                catch (Exception e)
                {
                    WriteLine(e.ToString(), nameof(io_proxy), nameof(ReadAllLines));
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        public static void WriteAllLines(string filename, IEnumerable<string> lines, string module_name="", string function_name="")
        {
            io_proxy.WriteLine($"{module_name}.{function_name} -> ( {filename} )", nameof(io_proxy), nameof(WriteAllLines));

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

        public static void AppendAllLines(string filename, IEnumerable<string> lines, string module_name="", string function_name="")
        {
            io_proxy.WriteLine($"{module_name}.{function_name} -> ( {filename} )", nameof(io_proxy), nameof(AppendAllLines));

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

        public static void AppendAllText(string filename, string text, string module_name="", string function_name="")
        {
            io_proxy.WriteLine($"{module_name}.{function_name} -> ( {filename} )", nameof(io_proxy), nameof(AppendAllText));

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

        public static void WriteAllText(string filename, string text, string module_name="", string function_name="")
        {
            io_proxy.WriteLine($"{module_name}.{function_name} -> ( {filename} )", nameof(io_proxy), nameof(WriteAllText));

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
