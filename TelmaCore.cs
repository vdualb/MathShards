#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

namespace MathShards.TelmaCore;

public enum AngleMeasureUnits { amuRadians = 0, amuDegrees = 1 };

public readonly struct PairReal : IEquatable<PairReal>
{
    public static readonly PairReal Zero = new PairReal(0, 0);
    public static readonly PairReal XAxis = new PairReal(1, 0);
    public static readonly PairReal YAxis = new PairReal(0, 1);

    public Real X { get; }
    public Real Y { get; }
    public Real[] AsArray() => new[] { X, Y };

    public PairReal(Real x, Real y)
    {
        X = x;
        Y = y;
    }

    public void Deconstruct(out Real x, out Real y)
        => (x, y) = (X, Y);
    /// <summary>
    ///  из полярных координат в декартовы
    /// </summary>
    public PairReal(Real a, Real b, AngleMeasureUnits measure) : this()
    {
        if (measure == AngleMeasureUnits.amuRadians)
        {
            X = a * (Real)Math.Cos(b);
            Y = a * (Real)Math.Sin(b);
        }
        else
        {
            Real c = b * (Real)Math.PI / 180;
            X = a * (Real)Math.Cos(c);
            Y = a * (Real)Math.Sin(c);
        }
    }
    public PairReal(ReadOnlySpan<Real> arr)
    {
#if DEBUG
        if (arr.Length != 2) throw new ArgumentException("Array size error");
#endif
        X = arr[0];
        Y = arr[1];
    }
    public Real this[int k]
    {
        get
        {
            return k switch
            {
                0 => X,
                1 => Y,
                _ => throw new Exception("get: Vector2D out of range"),
            };
        }
    }
    public static Real Distance(PairReal a, PairReal b) => (a - b).Norm;

    public static Real SqrDistance(PairReal a, PairReal b)
    {
        PairReal diff = a - b;
        return diff * diff;
    }
    public Real Distance(PairReal b) => (this - b).Norm;

    public Real SqrDistance(PairReal b) => SqrDistance(this, b);

    public Real Norm => (Real)Math.Sqrt(X * X + Y * Y);

    public PairReal Normalize() => this / Norm;

    public override string ToString() => $"Vec({X}, {Y})";

    public override bool Equals(object? obj) => obj is PairReal v && Equals(v);

    public override int GetHashCode() => HashCode.Combine(X, Y);

    public bool Equals(PairReal a) => a.X == X && a.Y == Y;
    public static bool TryParse(string line, out PairReal res)
    {
        Real x, y;
        var words = line.Split(new[] { ' ', '\t', ',', '>', '<', '(', ')' }, StringSplitOptions.RemoveEmptyEntries);
        if (words[0] == "Vec")
        {
            if (words.Length != 3 || !Real.TryParse(words[1], out x) || !Real.TryParse(words[2], out y))
            {
                res = Zero;
                return false;
            }
            else { res = new PairReal(x, y); return true; }
        }
        if (words.Length != 2 || !Real.TryParse(words[0], out x) || !Real.TryParse(words[1], out y))
        {
            res = Zero;
            return false;
        }
        else { res = new PairReal(x, y); return true; }
    }

    public static PairReal Parse(string line)
    {
        if (!TryParse(line, out PairReal res))
            throw new FormatException("Can't parse PairF64!");
        return res;
    }
    public PairReal Round(int digits) => new PairReal((Real)Math.Round(X, digits), (Real)Math.Round(Y, digits));

    public static PairReal Vec(Real x, Real y) => new PairReal(x, y);
    #region Static operators

    public static PairReal operator -(PairReal a) => new PairReal(-a.X, -a.Y);

    public static PairReal operator +(PairReal a, PairReal b) => new PairReal(a.X + b.X, a.Y + b.Y);

    public static PairReal operator -(PairReal a, PairReal b) => new PairReal(a.X - b.X, a.Y - b.Y);

    public static PairReal operator /(PairReal a, Real v) => new PairReal(a.X / v, a.Y / v);

    public static PairReal operator *(PairReal a, Real v) => new PairReal(a.X * v, a.Y * v);

    public static PairReal operator *(Real v, PairReal a) => new PairReal(v * a.X, v * a.Y);

    public static Real operator *(PairReal a, PairReal b) => a.X * b.X + a.Y * b.Y;

    public static bool operator ==(PairReal a, PairReal b) => a.X == b.X && a.Y == b.Y;

    public static bool operator !=(PairReal a, PairReal b) => a.X != b.X || a.Y != b.Y;

    public static PairReal Cross(PairReal v1) => new PairReal(v1.Y, -v1.X);

    public static Real Mixed(PairReal v1, PairReal v2) => v1.Y * v2.X - v1.X * v2.Y;

    public static PairReal Sum(PairReal a, PairReal b) => new PairReal(a.X + b.X, a.Y + b.Y);

    #endregion
    #region EqualityComparer

    private class EqualityComparer : IEqualityComparer<PairReal>
    {
        public int Digits { get; set; }

        public bool Equals(PairReal v1, PairReal v2)
        {
            return v1.Round(Digits) == v2.Round(Digits);
        }

        public int GetHashCode(PairReal obj)
        {
            return obj.Round(Digits).GetHashCode();
        }
    }

    public static IEqualityComparer<PairReal> CreateComparer(int digits = 7)
    {
        return new EqualityComparer { Digits = digits };
    }


    #endregion
}
