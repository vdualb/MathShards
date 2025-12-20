namespace MathShards.TelmaCore;

public enum AngleMeasureUnits { amuRadians = 0, amuDegrees = 1 };

public readonly struct PairF64 : IEquatable<PairF64>
{
    public static readonly PairF64 Zero = new PairF64(0, 0);
    public static readonly PairF64 XAxis = new PairF64(1, 0);
    public static readonly PairF64 YAxis = new PairF64(0, 1);

    public double X { get; }
    public double Y { get; }
    public double[] AsArray() => new[] { X, Y };

    public PairF64(double x, double y)
    {
        X = x;
        Y = y;
    }

    public void Deconstruct(out double x, out double y)
        => (x, y) = (X, Y);
    /// <summary>
    ///  из полярных координат в декартовы
    /// </summary>
    public PairF64(double a, double b, AngleMeasureUnits measure) : this()
    {
        if (measure == AngleMeasureUnits.amuRadians)
        {
            X = a * Math.Cos(b);
            Y = a * Math.Sin(b);
        }
        else
        {
            double c = b * Math.PI / 180;
            X = a * Math.Cos(c);
            Y = a * Math.Sin(c);
        }
    }
    public PairF64(ReadOnlySpan<double> arr)
    {
#if DEBUG
        if (arr.Length != 2) throw new ArgumentException("Array size error");
#endif
        X = arr[0];
        Y = arr[1];
    }
    public double this[int k]
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
    public static double Distance(PairF64 a, PairF64 b) => (a - b).Norm;

    public static double SqrDistance(PairF64 a, PairF64 b)
    {
        PairF64 diff = a - b;
        return diff * diff;
    }
    public double Distance(PairF64 b) => (this - b).Norm;

    public double SqrDistance(PairF64 b) => SqrDistance(this, b);

    public double Norm => Math.Sqrt(X * X + Y * Y);

    public PairF64 Normalize() => this / Norm;

    public override string ToString() => $"Vec({X}, {Y})";

    public override bool Equals(object? obj) => obj is PairF64 v && Equals(v);

    public override int GetHashCode() => HashCode.Combine(X, Y);

    public bool Equals(PairF64 a) => a.X == X && a.Y == Y;
    public static bool TryParse(string line, out PairF64 res)
    {
        double x, y;
        var words = line.Split(new[] { ' ', '\t', ',', '>', '<', '(', ')' }, StringSplitOptions.RemoveEmptyEntries);
        if (words[0] == "Vec")
        {
            if (words.Length != 3 || !double.TryParse(words[1], out x) || !double.TryParse(words[2], out y))
            {
                res = Zero;
                return false;
            }
            else { res = new PairF64(x, y); return true; }
        }
        if (words.Length != 2 || !double.TryParse(words[0], out x) || !double.TryParse(words[1], out y))
        {
            res = Zero;
            return false;
        }
        else { res = new PairF64(x, y); return true; }
    }

    public static PairF64 Parse(string line)
    {
        if (!TryParse(line, out PairF64 res))
            throw new FormatException("Can't parse PairF64!");
        return res;
    }
    public PairF64 Round(int digits) => new PairF64(Math.Round(X, digits), Math.Round(Y, digits));

    public static PairF64 Vec(double x, double y) => new PairF64(x, y);
    #region Static operators

    public static PairF64 operator -(PairF64 a) => new PairF64(-a.X, -a.Y);

    public static PairF64 operator +(PairF64 a, PairF64 b) => new PairF64(a.X + b.X, a.Y + b.Y);

    public static PairF64 operator -(PairF64 a, PairF64 b) => new PairF64(a.X - b.X, a.Y - b.Y);

    public static PairF64 operator /(PairF64 a, double v) => new PairF64(a.X / v, a.Y / v);

    public static PairF64 operator *(PairF64 a, double v) => new PairF64(a.X * v, a.Y * v);

    public static PairF64 operator *(double v, PairF64 a) => new PairF64(v * a.X, v * a.Y);

    public static double operator *(PairF64 a, PairF64 b) => a.X * b.X + a.Y * b.Y;

    public static bool operator ==(PairF64 a, PairF64 b) => a.X == b.X && a.Y == b.Y;

    public static bool operator !=(PairF64 a, PairF64 b) => a.X != b.X || a.Y != b.Y;

    public static PairF64 Cross(PairF64 v1) => new PairF64(v1.Y, -v1.X);

    public static double Mixed(PairF64 v1, PairF64 v2) => v1.Y * v2.X - v1.X * v2.Y;

    public static PairF64 Sum(PairF64 a, PairF64 b) => new PairF64(a.X + b.X, a.Y + b.Y);

    #endregion
    #region EqualityComparer

    private class EqualityComparer : IEqualityComparer<PairF64>
    {
        public int Digits { get; set; }

        public bool Equals(PairF64 v1, PairF64 v2)
        {
            return v1.Round(Digits) == v2.Round(Digits);
        }

        public int GetHashCode(PairF64 obj)
        {
            return obj.Round(Digits).GetHashCode();
        }
    }

    public static IEqualityComparer<PairF64> CreateComparer(int digits = 7)
    {
        return new EqualityComparer { Digits = digits };
    }


    #endregion
}
