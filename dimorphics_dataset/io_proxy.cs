using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class io_proxy
    {
        private static readonly object _console_lock = new object();

        internal static void WriteLine(string text = @"", string module_name = @"", string method_name = @"", bool use_lock = false)
        {
            void write(string text)
            {
                if (string.IsNullOrEmpty(text)) Console.WriteLine();
                else
                {
                    var pid = Environment.ProcessId;
                    var thread_id = Thread.CurrentThread.ManagedThreadId;
                    var task_id = Task.CurrentId ?? 0;

                    Console.WriteLine($@"{DateTime.Now:G} {pid:000000}.{thread_id:000000}.{task_id:000000} {module_name}.{method_name} -> {text}");
                }
            }

            if (!program.verbose) return;

            try
            {
                if (use_lock)
                {
                    lock (_console_lock)
                    {
                        write(text);
                    }
                }
                else
                {
                    write(text);
                }
            }
            catch (Exception)
            {

            }
        }

        internal static void CreateDirectory(string filename)//, string module_name, string method_name)
        {
            //program.WriteLine($@"{module_name}.{method_name} -> {nameof(program)}.{nameof(CreateDirectory)} ( {filename} )");

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

        internal static string[] ReadAllLines(string filename, string module_name="", string method_name="")
        {
            io_proxy.WriteLine($@"{module_name}.{method_name} -> ( {filename} )", nameof(io_proxy), nameof(ReadAllLines));

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

        internal static void WriteAllLines(string filename, IEnumerable<string> lines, string module_name="", string method_name="")
        {
            io_proxy.WriteLine($@"{module_name}.{method_name} -> ( {filename} )", nameof(io_proxy), nameof(WriteAllLines));

            CreateDirectory(filename);//, module_name, method_name);

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

        internal static void AppendAllLines(string filename, IEnumerable<string> lines, string module_name="", string method_name="")
        {
            io_proxy.WriteLine($@"{module_name}.{method_name} -> ( {filename} )", nameof(io_proxy), nameof(AppendAllLines));

            CreateDirectory(filename);//, module_name, method_name);

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

        internal static void AppendAllText(string filename, string text, string module_name="", string method_name="")
        {
            io_proxy.WriteLine($@"{module_name}.{method_name} -> ( {filename} )", nameof(io_proxy), nameof(AppendAllText));

            CreateDirectory(filename);//, module_name, method_name);

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

        internal static void WriteAllText(string filename, string text, string module_name="", string method_name="")
        {
            io_proxy.WriteLine($@"{module_name}.{method_name} -> ( {filename} )", nameof(io_proxy), nameof(WriteAllText));

            CreateDirectory(filename);//, module_name, method_name);

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
