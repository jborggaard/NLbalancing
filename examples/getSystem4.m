function [A,B,C,N,zInit,M] = getSystem4(n,m,p,epsilon)
%getSystem4  An input-output system based on discretization of the K-S equation.
%
%  Usage: [A,B,C,N,zInit] = getSystem4(n,m,p,epsilon)
%
%  n   is the order of the system (number of internal nodes)
%  m   is the number of equally spaced control inputs
%  p   is the number of observations
%
%  Produces a finite element model of the conservation form of the KS equations
%
%      z_t = -\epsilon z_xx - \epsilon^2 z_xxxx - 2\epsilon z z_x + 
%             \sum_{i=1}^m \chi_i(x) u_i(t)
%     y(t) = \int_0^1 \chi_i(x) z(x,t) dx
%
%  with boundary and initial conditions 
%
%      z(0,t)=z(1,t), z_x(0,t)=z_x(1,t), z(x,0) = z_0(x) = sin(4*pi*x)
%
%  These are discretized using piecewise Hermite cubic (H2) elements.
%
%  B represents (m) equispaced piecewise-constant source terms, and
%  C represents integral averages over (p) equispaced portions of the domain.
%
%  Default values:  n=32, m=1, epsilon=0.005890757777002
%                   (see Holmes, Lumley, and Berkooz)
%
%  A finite element model is generated, then placed into standard ODE form
%  (identity mass matrix) using a change of variables.  The final discretized
%  system has the form
%
%       z_t = Az + N kron(z,z) + Bu,   z(0) = zInit
%         y = Cz
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  References: Nonlinear balanced truncation: Part 1-Computing energy functions,
%              Kramer, Gugercin, and Borggaard.  arXiv
%
%              Nonlinear balanced truncation: Part 2-Model reduction on manifolds,
%              Kramer, Gugercin, and Borggaard.  arXiv
%
%              Nonlinear feedback control of PDEs with polynomial nonlinearities
%              IFAC-Online, 
%              Jeff Borggaard and Lizette Zietsman, 2021.
%
%  Part of the NLbalancing and QQR repositories.
%-------------------------------------------------------------------------------
%%

  if ( nargin<2 )
    n = 32;
    m = 1;
    p = 2;
  end
  
  if ( nargin<3 )
    epsilon = 1/(13.0291)^2;    % L = 13.0291 
  end
  L = 1/sqrt(epsilon);

  [x,eConn] = oned_mesh([0; 1],[1 2],n);
  
  [n_nodes   , ~       ] = size(x    );
  [n_elements, nel_dof ] = size(eConn);
  
  for n_nd=1:n_nodes-1
    ide(n_nd,1) = 2*n_nd-1;
    ide(n_nd,2) = 2*n_nd;
  end

  %  set the periodic conditions
  ide(n_nodes,1) = 1;
  ide(n_nodes,2) = 2;

  n_equations = 2*(n_nodes-1);

  %  Define an index function into N2
  tnm1 = 2*(n_nodes-1);
  idx2 = @(j,k) (j-1)*tnm1 + k;
  
  n_gauss = 5;    % number of quadrature points.
  [r,wt] = oned_quadrature(n_gauss);

  one    = ones(n_gauss,1);
  eps_g  = epsilon  *one;
  eps2_g = epsilon^2*one;

  II  = zeros(n_elements*2*nel_dof,1);
  JJ  = zeros(n_elements*2*nel_dof,1);
  AA  = zeros(n_elements*2*nel_dof,1);
  MM  = zeros(n_elements*2*nel_dof,1);
  B   = zeros(n_equations,m);
  C   = zeros(p,n_equations);
  z0  = zeros(n_equations,1);
  z0p = zeros(n_equations,1);

  IIn = zeros(4*n_elements*(2*nel_dof)^3,1);
  JJn = zeros(4*n_elements*(2*nel_dof)^3,1);
  NN  = zeros(4*n_elements*(2*nel_dof)^3,1);

  n_triplets = 0;
  n_tripletsn = 0;

  B0_loc = zeros(nel_dof,m);
  B1_loc = zeros(nel_dof,m);
  C0_loc = zeros(p,nel_dof);
  C1_loc = zeros(p,nel_dof);

  for n_el=1:n_elements
    % compute value of each test function and spatial derivatives
    % at the integration points (x_g - Gauss points, wt_g - Gauss weights)
    nodes_local = eConn(n_el,:);
    x_local     = x(nodes_local,:);
    [x_g,wt_g, phi0,phi1, p0_x,p1_x, p0_xx,p1_xx] = ...
                                               oned_shapeherm(x_local,r,wt);

    M00_loc =  oned_bilinear( one, phi0, phi0, wt_g );
    M01_loc =  oned_bilinear( one, phi0, phi1, wt_g );
    M10_loc =  M01_loc.';
    M11_loc =  oned_bilinear( one, phi1, phi1, wt_g );
    
    A00_loc =  oned_bilinear(  eps_g, p0_x, p0_x, wt_g ) ...
              -oned_bilinear( eps2_g, p0_xx, p0_xx, wt_g );
    A01_loc =  oned_bilinear(  eps_g, p0_x, p1_x, wt_g ) ...
              -oned_bilinear( eps2_g, p0_xx, p1_xx, wt_g );
    A10_loc =  A01_loc.';
    A11_loc =  oned_bilinear(  eps_g, p1_x, p1_x, wt_g ) ...
              -oned_bilinear( eps2_g, p1_xx, p1_xx, wt_g );
    
    for mm = 1:m
      b_loc = chi(x_g,mm,m);
      
      B0_loc(:,mm) = oned_f_int( b_loc, phi0, wt_g );
      B1_loc(:,mm) = oned_f_int( b_loc, phi1, wt_g );
    end

    for pp = 1:p
      c_loc = chi(x_g,pp,p);

      C0_loc(pp,:) = oned_f_int( c_loc, phi0, wt_g ).';
      C1_loc(pp,:) = oned_f_int( c_loc, phi1, wt_g ).';
    end

    [z_loc,~] = zZero(x_g,L);
    z0_loc    = oned_f_int( z_loc, phi0, wt_g );
    z1_loc    = oned_f_int( z_loc, phi1, wt_g );

    %---------------------------------------------------------------------------
    % Assemble contributions into the system matrices
    %---------------------------------------------------------------------------
    for n_t=1:nel_dof
      n_test0 = ide(nodes_local(n_t),1);
      n_test1 = ide(nodes_local(n_t),2);

      for n_u=1:nel_dof
        n_unk0 = ide(nodes_local(n_u),1);
        n_unk1 = ide(nodes_local(n_u),2);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk0;
        JJ(n_triplets) = n_test0;
        AA(n_triplets) = A00_loc(n_t,n_u);
        MM(n_triplets) = M00_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk0;
        JJ(n_triplets) = n_test1;
        AA(n_triplets) = A01_loc(n_t,n_u);
        MM(n_triplets) = M01_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk1;
        JJ(n_triplets) = n_test0;
        AA(n_triplets) = A10_loc(n_t,n_u);
        MM(n_triplets) = M10_loc(n_t,n_u);

        n_triplets = n_triplets + 1;
        II(n_triplets) = n_unk1;
        JJ(n_triplets) = n_test1;
        AA(n_triplets) = A11_loc(n_t,n_u);
        MM(n_triplets) = M11_loc(n_t,n_u);
      end
        
      for mm=1:m
        B(n_test0,mm) = B(n_test0,mm) + B0_loc(n_t,mm);
        B(n_test1,mm) = B(n_test1,mm) + B1_loc(n_t,mm);
      end
      
      for pp=1:p
        C(pp,n_test0) = C(pp,n_test0) + C0_loc(pp,n_t);
        C(pp,n_test1) = C(pp,n_test1) + C1_loc(pp,n_t);
      end
      
      z0(n_test0) = z0(n_test0) + z0_loc(n_t);
      z0(n_test1) = z0(n_test1) + z1_loc(n_t);
      
      %  assemble the 2zz_x term with symmetries
      for nj=1:nel_dof
        j0 = ide(nodes_local(nj),1);
        j1 = ide(nodes_local(nj),2);
        for nk=nj:nel_dof
          k0 = ide(nodes_local(nk),1);
          k1 = ide(nodes_local(nk),2);

          tmp00 = sum(p0_x(:,n_t).*phi0(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j0,k0);
          NN(n_tripletsn)  = tmp00;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k0,j0);
          NN(n_tripletsn)  = tmp00;

          tmp01 = sum(p0_x(:,n_t).*phi0(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j0,k1);
          NN(n_tripletsn)  = tmp01;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k1,j0);
          NN(n_tripletsn)  = tmp01;

          tmp10 = sum(p0_x(:,n_t).*phi1(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j1,k0);
          NN(n_tripletsn)  = tmp10;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k0,j1);
          NN(n_tripletsn)  = tmp10;

          tmp11 = sum(p0_x(:,n_t).*phi1(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(j1,k1);
          NN(n_tripletsn)  = tmp11;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test0;
          JJn(n_tripletsn) = idx2(k1,j1);
          NN(n_tripletsn)  = tmp11;

          %  for p1_x
          tmp00 = sum(p1_x(:,n_t).*phi0(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j0,k0);
          NN(n_tripletsn)  = tmp00;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k0,j0);
          NN(n_tripletsn)  = tmp00;

          tmp01 = sum(p1_x(:,n_t).*phi0(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j0,k1);
          NN(n_tripletsn)  = tmp01;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k1,j0);
          NN(n_tripletsn)  = tmp01;

          tmp10 = sum(p1_x(:,n_t).*phi1(:,nj).*phi0(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j1,k0);
          NN(n_tripletsn)  = tmp10;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k0,j1);
          NN(n_tripletsn)  = tmp10;

          tmp11 = sum(p1_x(:,n_t).*phi1(:,nj).*phi1(:,nk).*wt_g(:))/2;
          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(j1,k1);
          NN(n_tripletsn)  = tmp11;

          n_tripletsn = n_tripletsn+1;
          IIn(n_tripletsn) = n_test1;
          JJn(n_tripletsn) = idx2(k1,j1);
          NN(n_tripletsn)  = tmp11;

        end
      end
    end
  end

  II = II(1:n_triplets);
  JJ = JJ(1:n_triplets);
  A = sparse( II, JJ, AA(1:n_triplets), n_equations, n_equations );
  M = sparse( II, JJ, MM(1:n_triplets), n_equations, n_equations );
  
  N = sparse( IIn(1:n_tripletsn), JJn(1:n_tripletsn), NN(1:n_tripletsn), n_equations, n_equations^2 );
  N = epsilon*N;
  
  zInit = M\z0;

  %  We now have a system of the form
  %
  %   M \dot{z} = Az + Bu + N*kron(z,z)
  %
  %  Perform a change of variables to eliminate the positive definite
  %  mass matrix.  Let M^(1/2)z -> x
  %
  %    \dot{x} = Ax+Bu+N*kron(x,x),
  %
  %  This is required until we extend our formulation to handle "mass matrices"

  sqM = sqrtm(full(M));

  sqMinv = inv(sqM);

  A = sqMinv*A*sqMinv;                  %#ok
  B = sqMinv*B;                         %#ok
  C = C*sqMinv;                         %#ok
  N = kroneckerRight(sqMinv*N,sqMinv);  %#ok
  
  zInit = sqM*zInit;    % change of variable to remove mass matrix M

end

function [b_loc] = chi(x_local,mm,m)
%  The characteristic function over the interval ( (mm-1)/m, mm/m ).
  n = length(x_local);
  b_loc = zeros(n,1);
  
  for i=1:n
    if ( x_local(i)>(mm-1)/m && x_local(i)<mm/m )
      b_loc(i) = 1;
    else
      b_loc(i) = 0;
    end
  end
end

function [z0,z0p] = zZero(x,L)
%  Sets the initial function (FEM coefs are determined by projection)
  
  z0  = sin(4*pi*x)*L;
  z0p = 4*pi*cos(4*pi*x)*L;
  z0  = sin(pi*x)*L;
  z0p = pi*cos(pi*x)*L;
  
end

function M = oned_bilinear( kernel, phi, test, w_g )
%-----------------------------------------------------------------------
%  oned_bilinear.m - routine to compute \int{ kernel*phi*test }
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    M = oned_bilinear(kernel, phi, test, w_g)
%
%  Variables:     kernel
%                        Kernel function in the integral evaluated
%                        at the quadrature points
%
%                 phi
%                        matrix of element test functions evaluated
%                        at the quadrature points (dim: n_quadrature, n_dof)
%
%                 test
%                        matrix of test functions evaluated at the
%                        quadrature points (dim: n_quadrature, n_test)        
%
%                 w_g
%                        Column vector of quadrature weights
%-----------------------------------------------------------------------

  %  Vectorized version is more efficient (even for small vector lengths)
    M = test'*diag(kernel.*w_g)*phi;

end

function F = oned_f_int( Ff, test, w_g )
%-----------------------------------------------------------------------
%  oned_f_int.m - routine to compute \int{ f*test }
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    F = oned_f_int( Ff, test, w_g )
%
%  Variables:     Ff    
%                        Function values at the quadrature points
%
%                 test
%                        matrix of test functions evaluated at the
%                        quadrature points (dim: n_quadrature, n_dof)
%
%                 w_g
%                        Column vector of quadrature weights
%-----------------------------------------------------------------------

  %  Vectorized version is more efficient (even for small vector lengths)
  F = test'*(w_g.*Ff);

end

function [x,e_conn,index_u,index_c] = oned_mesh(xb,e_connb,rho)
%%----------------------------------------------------------------------
%  oned_mesh   - Generate a mesh with a prescribed density.
%                This routine returns elements of the same type as 
%                xb, e_connb (linear or quadratic)
%
%  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
%  Version: 1.0a
%
%  Usage:    [x,e_conn,index_u,index_c] = oned_mesh(xb,e_connb,rho)
%
%  Variables:     xb    
%                        nodal coordinates for a background mesh
%                 e_connb 
%                        connectivity for a background mesh
%                 rho     
%                        a mesh density function
%                        (assumed piecewise constant for now)
%
%                 x       
%                        node coordinates of adapted mesh
%                 e_conn  
%                        element connectivity of adapted mesh
%                 index_u 
%                        node numbers of unknowns
%                 index_c 
%                        node numbers of controls
%% ---------------------------------------------------------------------

  dim = size(e_connb,2);
  dim = dim - 1;

  rho = rho(:);  % make rho a column ;^)

  % make sure the number of elements is integral
  int_rho  = ( xb(e_connb(:,end),1)-xb(e_connb(:,1),1) )'*rho;
  new_elem = ceil(int_rho);
  rho      = new_elem*rho/int_rho;

  x_front  = xb(1,1);
  int_rho  = 0;
  bg_elem  = 0;

  for k=1:new_elem
    if (dim == 1)
      e_conn(k,:) = [k, k+1];
    elseif (dim == 2)
      e_conn(k,:) = [2*k-1, 2*k, 2*k+1];
    elseif (dim == 3)
      e_conn(k,:) = [3*k-2, 3*k-1, 3*k, 3*k+1];
    end 

    while (int_rho<1-sqrt(eps))
      bg_elem  = bg_elem + 1;
      eint_rho = ( xb(e_connb(bg_elem,end),1)-xb(e_connb(bg_elem,1),1) )*...
                 rho(bg_elem);
      int_rho  = int_rho + eint_rho;
    end

    % the new endpoint is in the current background element
    x_t     = max(x_front,xb(e_connb(bg_elem,1)));
    int_rho = int_rho-1;

    x_right = x_t + min(1,eint_rho-int_rho)/eint_rho*...
              ( xb(e_connb(bg_elem,end),1)-xb(e_connb(bg_elem,1),1) );
    x_nodes = linspace(x_front,x_right,dim+1);
    x(e_conn(k,:),1) = x_nodes';

    % advance the front
    x_front = x_right;
  end

  [n_nodes,~] = size(x);
  index_u = [2:n_nodes-1];
  index_c = [1, n_nodes];
end

function [r,w] = oned_quadrature(rule)
%-----------------------------------------------------------------------
%  oned_gauss.m - calculate Gauss integration points on (-1,1)
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [r,w] = oned_quadrature(rule)
%
%  Variables:     rule  
%                        Number of Gauss points:
%                 r
%                        Gauss points located between (-1,1)      
%                 w
%                        Gauss weights corresponding to r
%-----------------------------------------------------------------------

  r = zeros(rule,1);
  w = zeros(rule,1);

  if (rule == 1)       % up to order 1 polynomials exact
    r(1) = 0;
    w(1) = 2;
    
  elseif (rule == 2)   % up to order 3 polynomials exact
    r(1) =-1.0d0 / sqrt(3.0d0);
    r(2) =-r(1);
    w(1) = 1.0;
    w(2) = 1.0;
    
  elseif (rule == 3)  % up to order 5 polynomials exact
    r(1) =-sqrt(3.0d0/5.0d0);
    r(2) = 0.0;
    r(3) =-r(1);
    w(1) = 5.0d0 / 9.0d0;
    w(2) = 8.0d0 / 9.0d0;
    w(3) = w(1);
    
  elseif (rule == 4)  % up to order 7 polynomials exact
    r(1) =-sqrt((3.0d0+2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(2) =-sqrt((3.0d0-2.0*sqrt(6.0d0/5.0d0))/7.0d0);
    r(3) =-r(2);
    r(4) =-r(1);
    w(1) = 0.5d0 - 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(2) = 0.5d0 + 1.0d0 / ( 6.0d0 * sqrt(6.0d0/5.0d0) );
    w(3) = w(2);
    w(4) = w(1);
    
  elseif (rule == 5)  % up to order 9 polynomials exact
    r(1) =-sqrt(5.0d0+4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(2) =-sqrt(5.0d0-4.0d0*sqrt(5.0d0/14.0d0)) / 3.0d0;
    r(3) = 0.0d0;
    r(4) =-r(2);
    r(5) =-r(1);
    w(1) = 161.0d0/450.0d0-13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(2) = 161.0d0/450.0d0+13.0d0/(180.d0*sqrt(5.0d0/14.0d0));
    w(3) = 128.0d0/225.0d0;
    w(4) = w(2);
    w(5) = w(1);
    
  elseif (rule == 6)
    r(1) = -0.2386191860831969;
    r(2) = -0.6612093864662645;
    r(3) = -0.9324695142031521;
    r(4) = - r(1);
    r(5) = - r(2);
    r(6) = - r(3);
    w(1) = 0.4679139345726910;
    w(2) = 0.3607615730481386;
    w(3) = 0.1713244923791704;
    w(4) = w(1);
    w(5) = w(2);
    w(6) = w(3);
    
  elseif (rule == 7)
    r(1) = -0.9491079123427585;
    r(2) = -0.7415311855993945;
    r(3) = -0.4058451513773972;
    r(4) =  0.0000000000000000;
    r(5) = - r(3);
    r(6) = - r(2);
    r(7) = - r(1);
    w(1) = 0.1294849661688697;
    w(2) = 0.2797053914892766;
    w(3) = 0.3818300505051189;
    w(4) = 0.4179591836734694;
    w(5) = w(3);
    w(6) = w(2);
    w(7) = w(1);
    
  elseif (rule == 8)
    r(1) = -0.9602898564975363;
    r(2) = -0.7966664774136267;
    r(3) = -0.5255324099163290;
    r(4) = -0.1834346424956498;
    r(5) = - r(4);
    r(6) = - r(3);
    r(7) = - r(2);
    r(8) = - r(1);
    w(1) = 0.1012285362903763;
    w(2) = 0.2223810344533745;
    w(3) = 0.3137066458778873;
    w(4) = 0.3626837833783620;
    w(5) = w(4);
    w(6) = w(3);
    w(7) = w(2);
    w(8) = w(1);

  elseif (rule == 9)
    r(1) = -0.9681602395076261;
    r(2) = -0.8360311073266358;
    r(3) = -0.6133714327005904;
    r(4) = -0.3242534234038089;
    r(5) =  0.0000000000000000;
    r(6) = - r(4);
    r(7) = - r(3);
    r(8) = - r(2);
    r(9) = - r(1);
    w(1) = 0.0812743883615744;
    w(2) = 0.1806481606948574;
    w(3) = 0.2606106964029354;
    w(4) = 0.3123470770400029;
    w(5) = 0.3302393550012598;
    w(6) = w(4);
    w(7) = w(3);
    w(8) = w(2);
    w(9) = w(1);
  
  elseif (rule == 10)
    r( 1) = -0.9739065285171717;
    r( 2) = -0.8650633666889845;
    r( 3) = -0.6794095682990244;
    r( 4) = -0.4333953941292472;
    r( 5) = -0.1488743389816312;
    r( 6) = - r(5);
    r( 7) = - r(4);
    r( 8) = - r(3);
    r( 9) = - r(2);
    r(10) = - r(1);
    w( 1) = 0.0666713443086881;
    w( 2) = 0.1494513491505806;
    w( 3) = 0.2190863625159820;
    w( 4) = 0.2692667193099963;
    w( 5) = 0.2955242247147529;
    w( 6) = w(5);
    w( 7) = w(4);
    w( 8) = w(3);
    w( 9) = w(2);
    w(10) = w(1);

  elseif (rule == 11)
    r( 1) = -0.9782286581460570;
    r( 2) = -0.8870625997680953;
    r( 3) = -0.7301520055740494;
    r( 4) = -0.5190961292068118;
    r( 5) = -0.2695431559523450;
    r( 6) =  0.0000000000000000;
    r( 7) = - r(5);
    r( 8) = - r(4);
    r( 9) = - r(3);
    r(10) = - r(2);
    r(11) = - r(1);
    w( 1) = 0.0556685671161737;
    w( 2) = 0.1255803694649046;
    w( 3) = 0.1862902109277343;
    w( 4) = 0.2331937645919905;
    w( 5) = 0.2628045445102467;
    w( 6) = 0.2729250867779006;
    w( 7) = w(5);
    w( 8) = w(4);
    w( 9) = w(3);
    w(10) = w(2);
    w(11) = w(1);

  else
    error('Quadrature rule not supported')
    keyboard
  end

end

function [x_g,w_g,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = ...
         oned_shapeherm(x,r,w)
%-----------------------------------------------------------------------
%  oned_shapeherm.m - computes test functions and derivatives on a
%                     (C1+) Hermite element given element coordinates 
%                     and Gauss points.
%
%  Copyright (c) 2013, Jeff Borggaard, Virginia Tech
%  Version: 1.3
%
%  Usage:    [x_g,w_g,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = ...
%                                          oned_shapeherm(x_local,r,w)
%
%  Variables:     x_local
%                        Coordinates of the element nodes
%                 r
%                        Coordinates of Gauss points in (-1,1)
%                 w
%                        Gauss weights associated with r
%
%                 x_g
%                        Coordinates of Gauss points in the element
%                 w_g
%                        Gauss weights scaled by the element Jacobian
%                 phi*
%                        Value of element shape functions at x_g
%                 p*_x
%                        First spatial derivatives of phi
%                 p*_xx
%                        Second spatial derivatives of phi
%-----------------------------------------------------------------------

  % [n,t1] = size(x);   % n = 2, representing the endpoints of the interval
                        % t1 had better be 1, since this is ONED_shapeherm
  rule   = length(r);   % derive the order of the quadrature rule

  % Transform coordinates for linear elements
  len = ( x(2,1)-x(1,1) );
  c0 = len/2;
  c1 = ( x(2,1)+x(1,1) )/2;

  x_g = c0*r + c1;

  phi0(:,1) = ( 1-r ).^2 .* ( .25*r + .50);
  phi0(:,2) = ( 1+r ).^2 .* (-.25*r + .50);

  phi1(:,1) = ( 1-r ).^2 .* ( .125*r + .125)*len;
  phi1(:,2) = ( 1+r ).^2 .* ( .125*r - .125)*len;
    
  p0_r(:,1) = ( 1-r ) .* (-.75*r - .75);
  p0_r(:,2) = ( 1+r ) .* (-.75*r + .75);
    
  p1_r(:,1) = ( 1-r ) .* (-.375*r - .125)*len;
  p1_r(:,2) = ( 1+r ) .* ( .375*r - .125)*len;

  dxdr = c0;
  djac = dxdr;
  drdx = 1./djac;
    
    
  p0_x(:,1) = p0_r(:,1).*drdx;
  p0_x(:,2) = p0_r(:,2).*drdx;
    
  p1_x(:,1) = p1_r(:,1).*drdx;
  p1_x(:,2) = p1_r(:,2).*drdx;
    
  p0_rr(:,1) = 1.5*r;
  p0_rr(:,2) =-1.5*r;
    
  p1_rr(:,1) = (.75*r - .25)*len;
  p1_rr(:,2) = (.75*r + .25)*len;
    
  p0_xx(:,1) = p0_rr(:,1)*drdx^2;
  p0_xx(:,2) = p0_rr(:,2)*drdx^2;
    
  p1_xx(:,1) = p1_rr(:,1)*drdx^2;
  p1_xx(:,2) = p1_rr(:,2)*drdx^2;
    
  w_g = djac.*w;

end
