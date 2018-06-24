/* 
  uniform_smooth_noise64.hh
  
  Copyright (c) 2018 Ryan P. Nicholl, "Exaeta", <exaeta@protonmail.com>
  
    https://github.com/Exaeta/det_noise64
    
    uniform_smooth_noise64 is a N-dimensional uniform smooth noise 
    function, but I'm going to change the hash function into something 
    faster/more well known in a later version (and/or allow templating 
    the hash used). The next version or version with templated hash will
    be cross platform deterministic.
    
    It doesn't rely on any third party libraries except boost::multiprecision,
    which can be substituited for any library that provides a 128-bit 
    unsigned integer.
    
    Any donations would be appreciated. :)
  
  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#ifndef RPNX_G_NOISE_HH
#define RPNX_G_NOISE_HH



#include <array>
#include <cinttypes>
#include <random>
#include <boost/multiprecision/cpp_int.hpp> 

namespace rpnx
{
  
  using uint64_t = std::uint64_t;
  using uint128_t = boost::multiprecision::uint128_t;
  
  /*
    This function is probably not cryptographically secure.
    For insecure functions only.
  */
  template <std::size_t N>
  std::uint64_t det_point_noise64(std::array<std::uint64_t, N> inputs, std::size_t c = 4)
  {
    std::uint64_t output = 0xFAFAFAFAFAFAFAFAull;
    for (std::size_t k = 0; k < c; k++)
    {
      for (std::size_t i = 0; i < N; i++)
      {
        output = (output >> 17) | (output  << (64 - 17));
        output += inputs[i];
        output = (output >> ((k+inputs[i]) &0x0F)) | (output  << (64 - ((k+inputs[i]) &0x0F)));
        output ^= 7;
        output ^= ((output & (0xbeull << 31)) >> 21);
        output += 1;
        
      }
    }
    output = (output >> 17) | (output  << (64 - 17));
    
    return output;
  }
  
  
  
  template <std::size_t N>
  bool det_point_noise1(std::array<std::uint64_t, N> inputs)
  {
    
    std::uint64_t output = 0xFAFAFAFAFAFAFAFAull;
    for (std::size_t k = 0; k < 4; k++)
    {
      for (std::size_t i = 0; i < N; i++)
      {
        output = (output >> 17) | (output  << (64 - 17));
        output += inputs[i];      
        output ^= 7;
        output ^= ((output & (0xbeull << 31)) >> 21);
        output += k;
      }
    }
    output = (output >> 17) | (output  << (64 - 17));
    
    return 1 & ( (output >> 7) ^ (output >> 3) ^ (output >> 5) ^ 
    (output >> 2) ^ (output >> 11) ^ (output >> 17) ^ (output >> 13) );
  }

  
  /** A deterministic pseudo-random N-dimensional smooth noise function, 
      which, assuming the values of the input and hash function are 
      uniformly distributed, provides uniformly distributed output.
      
      The algorithm used by this function might be changed in a future
      version of this library.      
      
      @param C C is the template parameter that represents the number of
      bits used to interpolate between coordinates. A higher value of C
      will make the noise smooth over a larger range.
      
      @param N N is the number of elements in the array. Generally leave
      this to automatic deduction.
      
      @param inputs A std::array of N 64-bit integers
      
      @time O(N 2^N), O(1) as compiled
      
      @space O(2^N),  O(1) as compiled
  */
  template <std::size_t C = 32, std::size_t N = 1>
  std::uint64_t wave_noise64(std::array<std::uint64_t, N> inputs)
  {
    using u64 = std::uint64_t;
    using u128 = uint128_t;
    using std::uint8_t;
    
    u64 out=0;
    
    std::array<uint64_t, (1 << N)> field_corners;
    // 2 ^ N dimensions
    // gathers the values in each corner.
    
    std::mt19937 rg;
    
    std::uint8_t r = 0; // used later;
    
    std::array<std::uint64_t, N> field_inputs = inputs;
    // Each dimension has a set of inputs
      
    for (std::size_t k = 0; k < field_inputs.size(); k++)
    {
      field_inputs[k] >>= C;
    
      // only the first N-C bits are significant to the field
       
     
    }
      
    
    for (std::size_t i = 0; i < field_corners.size(); i++)
    {
   
      std::array<std::uint64_t, N> field_inputs2 = field_inputs;
      
      
      for (std::size_t k = 0; k < N; k++)
      {
        if (i & (1 << k)) field_inputs2[k]++;
        // We must adjust the field inputs depending on which "corners"
        // we are in.
      }
      
      field_corners[i] = det_point_noise64(field_inputs2, 8);
      // calculate the corner values
      
      r += field_corners[i] & 0xFF;
      // replace lost randomness with this value;
      
      
    }   
    
   // std::cout << "intial corner[0]=" << field_corners[0] << std::endl;
    
    std::array<std::uint64_t, N> sigvals;
    std::array<bool, N> bvals;
    for (std::size_t i = 0; i < N; i++)
    {
      sigvals[i] = inputs[i] & ((std::uint64_t(1) << C) - 1);
      bvals[i] = 1 & inputs[i] >> C;
    }
    // We need to know the dimensional significance of dimension in
    // order to interpolate the values correctly.
    
    
    // The following code interpolates the N-dimensional structure
    // in N log N time by collapsing each dimension one at a time
    // until the structure is 0-dimensional and has only 1 value
    
    std::size_t q = N;
    // q here is to avoid taking logs later
    
    
    out = 0;
    
    for (size_t z = 0; z < N; z++)
    {
      u64 a = det_point_noise64(std::array<u64, 1>{field_inputs[z]});
      u64 b = det_point_noise64(std::array<u64, 1>{field_inputs[z]+1});
      
      u64 v1 = sigvals[z];
      u64 v2 = (u64(1) << C) - v1;
      
      //std::swap(a, b);
      if (bvals[z]) 
      {
        std::swap(a, b);
        std::swap(v1, v2);
      }
      
      u64 m1 = static_cast<u64>( ((u128(b) * v1) + (u128(a) * v2)) >> C );
      
      out += 0;
      out += m1;
      
      
      
      
    }
    
     if (out & (std::uint64_t(1) << 63))
    {
      out = ((~out ^ 1) << 1) | 1;
    }
    else out = out << 1;
    
    /* The transformation above is supposed to redistribute the "out" 
    value such that values near 0 are "similar" so there aren't abrupt 
    transitions due to modular arithmetic wrapping from 1 to 0.
    */
    
    out ^= 1 & (r ^ (r >> 3) ^ (r >> 5));
    /*
      An unfortunate result of the above transformation is that although 
      it's uniform over the 64-bit space, the last bit only changes when
      the sign flips and thus effectively only 63 bits are adequately 
      random. To fix that some "random" bits are used.
    */
    
    return out;
    
  }
  
   
}

#endif
