using System;
using System.Windows.Forms;
using Fem2DFluid.Forms;

namespace Fem2DFluid
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }
}
