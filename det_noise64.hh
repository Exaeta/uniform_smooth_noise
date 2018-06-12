/* 
  det_noise64.hh
  
  Copyright (c) 2018 Ryan P. Nicholl, "Exaeta", <exaeta@protonmail.com>
  
    https://github.com/Exaeta/det_noise64
    
    The following code implements a deterministic N-dimensional smooth
    64-bit integer noise function. This function hasn't been tested extensively.
    
    Also the "det_point_noise64" function should probably be improved or
    replaced with something faster.
    
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
#include <boost/multiprecision/cpp_int.hpp> 

namespace ioncraft
{
  
  using uint64_t = std::uint64_t;
  using uint128_t = boost::multiprecision::uint128_t;
  
  /*
    This function is probably not cryptographically secure.
    For insecure functions only.
  */
  template <std::size_t N>
  std::uint64_t det_point_noise64(std::array<std::uint64_t, N> inputs, std::size_t c = 9)
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
  
  template <std::size_t C = 32, std::size_t N = 1>
  std::uint64_t det_field_noise64(std::array<std::uint64_t, N> inputs, std::size_t c = 9)
  {
    
    std::array<std::uint64_t, (1 << N)> field_corners;
    // 2 ^ N dimensions
    // gathers the values in each corner.
    
    std::uint8_t r = 0; // used later;
    
    
    for (std::size_t i = 0; i < field_corners.size(); i++)
    {
      std::array<std::uint64_t, N> field_inputs = inputs;
      // Each dimension has a set of inputs
      
      
      for (std::size_t k = 0; k < field_inputs.size(); k++)
      {
        field_inputs[k] >>= C;
        // only the first N-C bits are significant to the field
      }
      
      
      for (std::size_t k = 0; k < N; k++)
      {
        if (i & (1 << k)) field_inputs[k]++;
        // We must adjust the field inputs depending on which "corners"
        // we are in.
      }
      
      field_corners[i] = det_point_noise64(field_inputs, c);
      // calculate the corner values
      
      r += field_corners[i] & 0xFF;
      // extra randomness to replace randomness lost later
    }   
    
    std::array<std::uint64_t, N> sigvals;
    for (std::size_t i = 0; i < N; i++)
    {
      sigvals[i] = inputs[i] & ((std::uint64_t(1) << C) - 1);
    }
    // We need to know the dimensional significance of dimension in
    // order to interpolate the values correctly.
    
    
    // The following code interpolates the N-dimensional structure
    // in N log N time by collapsing each dimension one at a time
    // until the structure is 0-dimensional and has only 1 value
    
    std::size_t q = N;
    // q here is to avoid taking logs later
    
    for (std::size_t m = field_corners.size() >> 1
     /* The initial value of m here is half the number of corners, 
        because the collapsing inner
        loop assumes that i|m is valid.
        The number of corners will always be a power of two, so this works.
     */; 
     m != 0;
     m >>= 1
     /*
      When a round finishes, we can collapse the next dimension, which
      is represented by a power of two.
     */)
    {
      q--;
      r ++;
      // reduce by one before the loop (so sigvals[q] is valid because
      // q < N
      
      for (std::size_t i = 0; i < m; i++)
      {
        /* There are two corners involved in each inner loop of the
           collapse, e.g. (0, x, y, z) and (1, x, y, z) where x, y, and
           z are 0 or 1 and i = 0b[wxyz] 
           so (w << 3) | (x << 2) | (y << 1) | (z << 0) etc. but for
           N-dimensions.
           */
           
        std::uint64_t a = field_corners[i];
        // a is the first corner
        std::uint64_t b = field_corners[i | m];
        // and b is the second
        
        std::uint64_t v = sigvals[q];
        // a is the "low" value is a low value for v should indicate a high contribution for a
        // b is the high value so a high value of v should be a high contribution from b
        uint128_t o = uint128_t(a) * 2 * ((uint128_t(1) << (C)) - v) + 2 * v * uint128_t(b);
        
        r ^= static_cast<std::uint8_t>((field_corners[i] & 0xFF00) >> 8);
        // update: use boost's 128-bit number to get around precision loss.
        // It's unfortunate that precision is lost here... maybe there is a better way to do this
        
        field_corners[i] = static_cast<std::uint64_t> ((o >> (C-2)) & ((uint128_t(1) << 64)-1));
        field_corners[i] ^= 3 & ( (r >> 3) ^ (r >> 7));
      }
    }
    
    std::uint64_t out = field_corners[0];
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
      random. To fix that some random bits are used.
    */
    
    
    return out;
    
    
  }
}

#endif
