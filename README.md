# RPNX Uniform Smooth Noise

```
template <size_t C = 32, size_t N = 1>
std::uint64_t rpnx::uniform_smooth_noise64(std::array<std::uint64_t, N> inputs);
```

Higher C-value gives more bits of interpolation. (slower change in noise value)

## Example:

```
std::array<uint64_t, 3> xyz = {pos.x, pos.y, pos.z};
uint64_t pseudorandom_value = rpnx::uniform_smooth_noise<16>(xyz);
```
