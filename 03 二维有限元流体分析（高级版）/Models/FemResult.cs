using System;

namespace Fem2DFluid.Models
{
    [Serializable]
    public class FemResult
    {
        private bool _converged;
        private int _iterations;
        private long _assembleMs;
        private long _solveMs;
        private long _postMs;
        private double _maxPhi;
        private double _minPhi;
        private double _maxV;
        private double _minV;
        private string _message;

        public FemResult()
        {
            _converged = false;
            _iterations = 0;
            _assembleMs = 0;
            _solveMs = 0;
            _postMs = 0;
            _maxPhi = 0.0;
            _minPhi = 0.0;
            _maxV = 0.0;
            _minV = 0.0;
            _message = string.Empty;
        }

        public bool Converged
        {
            get { return _converged; }
            set { _converged = value; }
        }

        public int Iterations
        {
            get { return _iterations; }
            set { _iterations = value; }
        }

        public long AssembleMs
        {
            get { return _assembleMs; }
            set { _assembleMs = value; }
        }

        public long SolveMs
        {
            get { return _solveMs; }
            set { _solveMs = value; }
        }

        public long PostMs
        {
            get { return _postMs; }
            set { _postMs = value; }
        }

        public double MaxPhi
        {
            get { return _maxPhi; }
            set { _maxPhi = value; }
        }

        public double MinPhi
        {
            get { return _minPhi; }
            set { _minPhi = value; }
        }

        public double MaxV
        {
            get { return _maxV; }
            set { _maxV = value; }
        }

        public double MinV
        {
            get { return _minV; }
            set { _minV = value; }
        }

        public string Message
        {
            get { return _message; }
            set { _message = value; }
        }

        public long TotalMs
        {
            get { return _assembleMs + _solveMs + _postMs; }
        }

        public override string ToString()
        {
            return string.Format(
                "Result: converged={0}, iter={1}, asm={2}ms, sol={3}ms, post={4}ms, {5}",
                _converged, _iterations, _assembleMs, _solveMs, _postMs, _message ?? "");
        }
    }
}
