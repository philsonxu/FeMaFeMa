using System;
using System.Windows.Forms;
using Fem2DFluid.Forms;

namespace Fem2DFluid
{
    internal static class Program
    {
        [STAThread]
        private static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }
}
