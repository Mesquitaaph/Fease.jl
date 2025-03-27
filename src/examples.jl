function Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
  return (alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, a = a, b = b, u = u, u_x = u_x, f = f)
end

function example1()
  alpha = 1.0; beta = 1.0; gamma = 0.0; sigma(x) = 0; a = 0.0; b = 1.0;
  u(x) = x + (exp(-x) - exp(x))/(exp(1.0) - exp(-1));
  u_x(x) = 1 - (exp(-x) + exp(x)) *(1/(exp(1.0) - exp(-1)));
  f(x) = x

  return Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
end

function example2()
  alpha = pi; beta = exp(1.0); gamma = 0.0; sigma(x) = 0; a = 0.0; b = 1.0;
  u(x) = sin(pi*x);
  u_x(x) = cos(pi*x)*pi; u_xx(x) = -sin(pi*x)*pi^2;
  f(x) = -alpha*u_xx(x) + beta*u(x)

  return Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
end

function example3()
  alpha = 1.0; beta = 1.0; gamma = 1.0; sigma(x) = 0; a = 0.0; b = 1.0;
  u(x) = 4*(x - 1/2)^2 - 1;
  u_x(x) = 8*(x - 1/2); u_xx(x) = 8;
  f(x) = -alpha*u_xx(x) + gamma*u_x(x) + beta*u(x)

  return Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
end

function example4()
  alpha = 1.0; beta = 0.0; gamma = 0.0; sigma(x) = 1.0; a = 0.0; b = 1.0;
  u(x) = x + (exp(-x) - exp(x))/(exp(1.0) - exp(-1));
  u_x(x) = 1 - (exp(-x) + exp(x)) * (1/(exp(1.0) - exp(-1)));
  u_xx(x) = (1/(exp(1.0) - exp(-1))) * (-exp(x) + exp(-x));
  f(x) = -alpha*u_xx(x) + sigma(x)*u(x)

  return Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
end

function example5()
  alpha = 1.0; beta = 0.0; gamma = 0.0; sigma(x) = sin(x); a = 0.0; b = 1.0;
  u(x) = x + (exp(-x) - exp(x))/(exp(1.0) - exp(-1));
  u_x(x) = 1 - (exp(-x) + exp(x)) *(1/(exp(1.0) - exp(-1)));
  u_xx(x) = (1/(exp(1.0) - exp(-1))) * (-exp(x) + exp(-x));
  f(x) = -alpha*u_xx(x) + sigma(x)*u(x)

  return Example(alpha, beta, gamma, sigma, a, b, u, u_x, f)
end

function examples(case)
  if case == 1
      return example1()
  elseif case == 2
      return example2()
  elseif case == 3
      return example3()
  elseif case == 4
      return example4()
  elseif case == 5
      return example5()
  end
end