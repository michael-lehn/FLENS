namespace flens {

#ifdef AUX_ROUND
template <typename T>
T
round(T x)
{
    if (x>=T(0)) {
        T y = floor(x);
        if (x-y>=T(0.5)) {
            y += T(1);
        }
        return y;
    }
    T y = ceil(x);
    if (y-x>=T(0.5)) {
        y -= T(1.0);
    }
    return y;
}
#endif // AUX_ROUND

} // namespace flens
