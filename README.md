# det_noise64
Deterministic Smooth Noise Functions



Pass the following function a std::array of 64 bit integers.

```
  template <std::size_t C = 32, std::size_t N = 1>
  std::uint64_t det_field_noise64(std::array<std::uint64_t, N> inputs, std::size_t c = 9)
```

