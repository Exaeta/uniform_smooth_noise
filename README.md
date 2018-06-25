# RPNX Uniform Smooth Noise

Check out the examples:
https://github.com/Exaeta/uniform_smooth_noise/wiki/Examples


```
template <size_t C = 32, size_t N = 1>
std::uint64_t rpnx::wave_noise64(std::array<std::uint64_t, N> inputs);
```

Higher C-value gives more bits of interpolation. (slower change in noise value)

The wave noise function is smooth and uniform but has some directional artifacts and is (slightly) predictable. It's intended to be used as a low level function for building other noise functions.

## Example:

```
std::array<uint64_t, 3> xyz = {pos.x, pos.y, pos.z};
uint64_t pseudorandom_value = rpnx::wave_noise<16>(xyz);
```

## Help/Issues

Issues are tracked on github.
https://github.com/Exaeta/uniform_smooth_noise/issues
You can also get help using the library on irc.freenode.net, channel #rpnx, which is the preferred method unless you find a bug.
