using System;
using System.Collections.Generic;
using System.Text;

namespace Fem2DFluid.Models
{
    [Serializable]
    public class FemModel
    {
        private List<Node> _nodes;
        private List<Element> _elements;
        private List<Material> _materials;
        private HashSet<Edge> _boundaryEdges;
        private List<BoundaryEdge> _boundaryEdgeCache;
        private double _xmin;
        private double _xmax;
        private double _ymin;
        private double _ymax;
        private string _name;
        private bool _isNavierStokes;
        private double _dt;
        private int _timeStep;

        public FemModel()
        {
            _nodes = new List<Node>();
            _elements = new List<Element>();
            _materials = new List<Material>();
            _boundaryEdges = new HashSet<Edge>();
            _boundaryEdgeCache = new List<BoundaryEdge>();
            _xmin = 0.0;
            _xmax = 1.0;
            _ymin = 0.0;
            _ymax = 1.0;
            _name = "Model";
            _isNavierStokes = false;
            _dt = 0.01;
            _timeStep = 0;
            Material def = new Material(0, "Water");
            _materials.Add(def);
        }

        public List<Node> Nodes
        {
            get { return _nodes; }
            set { _nodes = value; }
        }

        public List<Element> Elements
        {
            get { return _elements; }
            set { _elements = value; }
        }

        public List<Material> Materials
        {
            get { return _materials; }
            set { _materials = value; }
        }

        public HashSet<Edge> BoundaryEdges
        {
            get { return _boundaryEdges; }
            set { _boundaryEdges = value; }
        }

        public List<BoundaryEdge> BoundaryEdgeCache
        {
            get { return _boundaryEdgeCache; }
            set { _boundaryEdgeCache = value; }
        }

        public double Xmin
        {
            get { return _xmin; }
            set { _xmin = value; }
        }

        public double Xmax
        {
            get { return _xmax; }
            set { _xmax = value; }
        }

        public double Ymin
        {
            get { return _ymin; }
            set { _ymin = value; }
        }

        public double Ymax
        {
            get { return _ymax; }
            set { _ymax = value; }
        }

        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }

        public bool IsNavierStokes
        {
            get { return _isNavierStokes; }
            set { _isNavierStokes = value; }
        }

        public double Dt
        {
            get { return _dt; }
            set { _dt = value; }
        }

        public int TimeStep
        {
            get { return _timeStep; }
            set { _timeStep = value; }
        }

        public Node GetNode(int id)
        {
            for (int i = 0; i < _nodes.Count; i++)
            {
                if (_nodes[i].Id == id)
                {
                    return _nodes[i];
                }
            }
            return null;
        }

        public Element GetElement(int id)
        {
            for (int i = 0; i < _elements.Count; i++)
            {
                if (_elements[i].Id == id)
                {
                    return _elements[i];
                }
            }
            return null;
        }

        public void ComputeBounds()
        {
            if (_nodes.Count == 0)
            {
                _xmin = 0.0; _xmax = 1.0; _ymin = 0.0; _ymax = 1.0;
                return;
            }
            _xmin = double.MaxValue;
            _xmax = -double.MaxValue;
            _ymin = double.MaxValue;
            _ymax = -double.MaxValue;
            for (int i = 0; i < _nodes.Count; i++)
            {
                Node n = _nodes[i];
                if (n.X < _xmin) _xmin = n.X;
                if (n.X > _xmax) _xmax = n.X;
                if (n.Y < _ymin) _ymin = n.Y;
                if (n.Y > _ymax) _ymax = n.Y;
            }
        }

        public void ResetResults()
        {
            for (int i = 0; i < _nodes.Count; i++)
            {
                Node n = _nodes[i];
                n.Phi = 0.0;
                n.U = 0.0;
                n.V = 0.0;
                n.Vx = 0.0;
                n.Vy = 0.0;
                n.P = 0.0;
            }
            for (int i = 0; i < _elements.Count; i++)
            {
                Element e = _elements[i];
                e.Vx = 0.0;
                e.Vy = 0.0;
                e.Vmag = 0.0;
                e.P = 0.0;
            }
            _timeStep = 0;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat("FemModel '{0}': {1} nodes, {2} elements, {3} boundary edges",
                _name, _nodes.Count, _elements.Count, _boundaryEdgeCache.Count);
            return sb.ToString();
        }
    }
}
