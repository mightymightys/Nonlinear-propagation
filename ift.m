function result = ift(data)
  %result = fftshift(ifft(fftshift(data)));   %M.Joffre's version
  result = fftshift(ifft(ifftshift(data)));
