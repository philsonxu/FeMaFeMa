// Program.cs - 程序入口
using System;
using System.Windows.Forms;

namespace FEM2D
{
    internal static class Program
    {
        [STAThread]
        private static void Main()
        {
            ApplicationConfiguration.Initialize();
            Application.Run(new MainForm());
        }
    }
}
