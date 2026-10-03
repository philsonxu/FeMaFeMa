using System;

namespace Fem2DFluid.Models
{
    [Serializable]
    public class Material
    {
        private int _id;
        private string _name;
        private double _k;
        private double _nu;
        private double _density;
        private double _viscosity;

        public Material()
        {
            _id = 0;
            _name = "Default";
            _k = 1.0;
            _nu = 0.3;
            _density = 1000.0;
            _viscosity = 1.0e-3;
        }

        public Material(int id, string name)
            : this()
        {
            _id = id;
            _name = name;
        }

        public int Id
        {
            get { return _id; }
            set { _id = value; }
        }

        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }

        /// <summary>Permeability / hydraulic conductivity for seepage / potential flow</summary>
        public double K
        {
            get { return _k; }
            set { _k = value; }
        }

        /// <summary>Poisson ratio (reserved for poroelastic coupling)</summary>
        public double Nu
        {
            get { return _nu; }
            set { _nu = value; }
        }

        /// <summary>Density rho (kg/m^3)</summary>
        public double Density
        {
            get { return _density; }
            set { _density = value; }
        }

        /// <summary>Dynamic viscosity mu (Pa·s)</summary>
        public double Viscosity
        {
            get { return _viscosity; }
            set { _viscosity = value; }
        }

        /// <summary>Kinematic viscosity nu = mu / rho</summary>
        public double KinematicViscosity
        {
            get
            {
                if (Math.Abs(_density) < 1.0e-15)
                {
                    return _viscosity;
                }
                return _viscosity / _density;
            }
        }

        public override string ToString()
        {
            return string.Format("Material({0},{1})", _id, _name ?? "");
        }
    }
}
