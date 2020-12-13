using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

namespace dimorphics_dataset
{
    internal static class io_proxy
    {
        public const string module_name = nameof(io_proxy);
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
                //io_proxy.log_exception(e, /*program.string_debug*/($@""), module_name, method_name);
            }
        }

        internal static void log_exception(Exception e, string msg, string caller_module_name, string caller_method_name)
        {
            do
            {
                io_proxy.WriteLine($@"Error: ""{msg}"" ""{e.GetType()}"" ""{e.Source}"" ""{e.Message}"" ""{e.StackTrace}""", caller_module_name, caller_method_name);

#if DEBUG
                if (e.InnerException == null || e == e.InnerException)
                {
                    throw e;
                }
#endif
                e = e != e?.InnerException ? e?.InnerException : null;
            } while (e != null);
        }

        internal static void CreateDirectory(string filename)//, string module_name, string method_name)
        {
            //program.WriteLine(/*program.string_debug*/($@"{module_name}.{method_name} -> {nameof(program)}.{nameof(CreateDirectory)} ( {filename} )");

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

        internal static string[] ReadAllLines(string filename, string caller_module_name="", string caller_method_name="")
        {
            
            const string method_name = nameof(ReadAllLines);

            io_proxy.WriteLine(/*program.string_debug*/($@"{caller_module_name}.{caller_method_name} -> ( {filename} )"), module_name, method_name);

            while (true)
            {
                try
                {
                    var ret = File.ReadAllLines(filename);

                    return ret;
                }
                catch (Exception e)
                {
                    io_proxy.log_exception(e, /*program.string_debug*/($@"{caller_module_name}.{caller_method_name}"), module_name, method_name);
                    
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        internal static void WriteAllLines(string filename, IEnumerable<string> lines, string caller_module_name="", string caller_method_name="")
        {
            
            const string method_name = nameof(WriteAllLines);

            io_proxy.WriteLine(/*program.string_debug*/($@"{caller_module_name}.{caller_method_name} -> ( {filename} )"), module_name, method_name);

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
                    io_proxy.log_exception(e, /*program.string_debug*/($@"{caller_module_name}.{caller_method_name}"), module_name, method_name);
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        internal static void AppendAllLines(string filename, IEnumerable<string> lines, string caller_module_name = "", string caller_method_name = "")
        {
            
            const string method_name = nameof(AppendAllLines);

            io_proxy.WriteLine(/*program.string_debug*/($@"{caller_module_name}.{caller_method_name} -> ( {filename} )"), module_name, method_name);

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
                    io_proxy.log_exception(e, /*program.string_debug*/($@"{caller_module_name}.{caller_method_name}"), module_name, method_name);
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        internal static void AppendAllText(string filename, string text, string caller_module_name = "", string caller_method_name = "")
        {
            
            const string method_name = nameof(AppendAllText);

            io_proxy.WriteLine(/*program.string_debug*/($@"{caller_module_name}.{caller_method_name} -> ( {filename} )"), module_name, method_name);

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
                    io_proxy.log_exception(e, /*program.string_debug*/($@"{caller_module_name}.{caller_method_name}"), module_name, method_name);
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }

        internal static void WriteAllText(string filename, string text, string caller_module_name = "", string caller_method_name = "")
        {
            
            const string method_name = nameof(WriteAllText);


            io_proxy.WriteLine(/*program.string_debug*/($@"{caller_module_name}.{caller_method_name} -> ( {filename} )"), module_name, method_name);

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
                    io_proxy.log_exception(e, /*program.string_debug*/($@"{caller_module_name}.{caller_method_name}"), module_name, method_name);
                    Task.Delay(new TimeSpan(0, 0, 0, 10)).Wait();
                }
            }
        }
    }
}
