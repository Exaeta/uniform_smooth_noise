# det_noise64
Uniform Smooth Noise Function

```
template <size_t C = 32, size_t N = 1>
std::uint64_t rpnx::uniform_smooth_noise64(std::array<std::uint64_t, N> inputs);
```

Higher C-value gives more bits of interpolation. (slower change in noise value)
