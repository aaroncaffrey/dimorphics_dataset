using System;
using System.Diagnostics;
using System.Runtime;
using System.Runtime.Loader;
using System.Threading;

namespace dimorphics_dataset
{
    internal class init
    {
        const string module_name = nameof(init);

        internal static void set_process_priority()
        {
            const string method_name = nameof(set_process_priority);
            try { Process.GetCurrentProcess().PriorityBoostEnabled = false; } catch (Exception e) { io_proxy.log_exception(e, $@"", module_name, method_name); }
            try { Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.Idle; } catch (Exception e) { io_proxy.log_exception(e, $@"", module_name, method_name); }
        }

        internal static void set_thread_counts()
        {
            const string method_name = nameof(set_thread_counts);

            ThreadPool.SetMinThreads(Environment.ProcessorCount * 10, Environment.ProcessorCount * 10);
            ThreadPool.SetMaxThreads(Environment.ProcessorCount * 100, Environment.ProcessorCount * 100);
        }

        internal static void close_notifications(CancellationTokenSource cts = null)
        {
            const string method_name = nameof(close_notifications);

            Console.CancelKeyPress += (sender, eventArgs) =>
            {
                io_proxy.WriteLine($@"Console.CancelKeyPress", module_name, method_name);
                cts?.Cancel();
            };
            AssemblyLoadContext.Default.Unloading += context =>
            {
                io_proxy.WriteLine($@"AssemblyLoadContext.Default.Unloading", module_name, method_name);
                cts?.Cancel();
            };
            AppDomain.CurrentDomain.ProcessExit += (sender, eventArgs) =>
            {
                io_proxy.WriteLine($@"AppDomain.CurrentDomain.ProcessExit", module_name, method_name);
                cts?.Cancel();
            };
        }

        internal static void check_x64()
        {
            const string method_name = nameof(check_x64);

            var is_x64 = IntPtr.Size == 8;

            if (!is_x64) { throw new Exception($"{module_name}.{method_name}: Must run in 64bit mode: {nameof(is_x64)}={is_x64}."); }
        }

        internal static void set_gc_mode()
        {
            const string method_name = nameof(set_gc_mode);

            GCSettings.LatencyMode = GCLatencyMode.Batch;
            //GCSettings.LatencyMode = GCLatencyMode.SustainedLowLatency;
        }
    }
}