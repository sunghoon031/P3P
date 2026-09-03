#pragma once
#include <algorithm>
#include <array>
#include <cmath>


namespace p3p {


template <class T, int refinement_iterations = 5>
class P3P {
private:
static double solve_cubic_one_real_root(double a, double b, double c) {
    // Returns one real root of the monic cubic equation
    //     x^3 + a*x^2 + b*x + c = 0.
    // Substituting x = y - a/3 eliminates the quadratic term and gives
    // the depressed cubic
    //     y^3 - 3*q*y + 2*r = 0,
    // where q and r are defined below.
    double q = (a*a - 3.0*b) / 9.0;
    double r = (a*(2.0*a*a - 9.0*b) + 27.0*c) / 54.0;

    // The sign of r^2 - q^3 determines the root structure:
    //   r^2 < q^3:  three distinct real roots
    //   r^2 = q^3:  repeated real roots
    //   r^2 > q^3:  one real root and two complex-conjugate roots
    const double r2 = r*r;
    const double q3 = q*q*q;

    double x;
    if (r2 < q3) {
        // When all three roots are real.
        // y = -2*sqrt(q)*cos(theta/3), where cos(theta) = r / q^(3/2).
        const double q_sqrt = std::sqrt(q);
        const double cos_theta = std::clamp(r / (q*q_sqrt), -1.0, 1.0);
        const double y = -2.0*q_sqrt*std::cos(std::acos(cos_theta) / 3.0);

        // Undo the substitution x = y - a/3.
        x = y - a/3.0;
    }
    else {
        // Cardano's formula for the case with one real root, including
        // the limiting case where some roots are repeated:
        //     A = -sign(r) * cbrt(|r| + sqrt(r^2 - q^3))
        //     B = q/A
        //     y = A + B.
        double A = -std::cbrt(std::abs(r)+std::sqrt(r2-q3));
        if (r < 0.0) A = -A;

        // If A is zero, then q is also zero in the exact degenerate case.
        // Avoid division by zero and take B = 0.
        const double B = (A == 0.0) ? 0.0 : q/A;

        // Undo the substitution x = y - a/3.
        x = A + B - a/3.0;
    }
    return x;
}

static void solve_quartic_FerrariLagrange(double a, double b, double c, double d,
                                double roots[4], int& n_roots) {
    // Finds the real roots of the monic quartic polynomial
    //     x^4 + a*x^3 + b*x^2 + c*x + d = 0
    // using the Ferrari–Lagrange factorization:
    //     (x^2 + p1*x + q1)(x^2 + p2*x + q2) = 0.
    // Expanding the product gives the coefficient constraints
    //     p1 + p2       = a,
    //     p1*p2 + q1+q2 = b,
    //     p1*q2 + p2*q1 = c,
    //     q1*q2         = d.
    // Define
    //     y = q1 + q2.
    // Eliminating p1, p2, q1, and q2 produces the resolvent cubic
    //     y^3 - b*y^2 + (a*c - 4*d)*y
    //         + (4*b*d - a^2*d - c^2) = 0.
    // A real solution y lets us reconstruct the two quadratic factors.

    const double y = solve_cubic_one_real_root(-b, a*c - 4.0*d, -a*a*d - c*c + 4.0*b*d);

    // Since q1 + q2 = y and q1*q2 = d, q1 and q2 are the
    // solutions of
    //     q^2 - y*q + d = 0.
    // Its discriminant is y^2 - 4*d.
    const double D = y*y - 4.0*d;
    if (D < 0.0) return;
    const double sqrtD = std::sqrt(D);
    if (sqrtD == 0.0) return;

    // Recover q1 and q2 using quadratic formula:
    const double q1 = (y + sqrtD)*0.5;
    const double q2 = (y - sqrtD)*0.5;

    // Recover p1 and p2 using
    //     p1 + p2       = a,
    //     p1*q2 + p2*q1 = c.
    // Because q1 - q2 = sqrt_q_discriminant, solving this
    // two-equation linear system gives the expressions below.
    const double p1 = (a*q1 - c)/sqrtD;
    const double p2 = (c - a*q2)/sqrtD;

    // Solve the first quadratic factor: x^2 + p1*x + q1 = 0.
    const double D1 = p1*p1 - 4.0*q1;
    if (D1 >= 0.0) {
        const double s = std::sqrt(D1);
        roots[n_roots++] = (s-p1)*0.5;
        roots[n_roots++] = (-s-p1)*0.5;
    }

    // Solve the second quadratic factor: x^2 + p2*x + q2 = 0.
    const double D2 = p2*p2 - 4.0*q2;
    if (D2 >= 0.0) {
        const double s = std::sqrt(D2);
        roots[n_roots++] = (s-p2)*0.5;
        roots[n_roots++] = (-s-p2)*0.5;
    }
}

static void solve_quartic_ClassicalFerrari(double a, double b, double c, double d,
                                           double roots[4], int& n_roots) {
    // Finds the real roots of the monic quartic equation
    //     x^4 + a*x^3 + b*x^2 + c*x + d = 0
    // using the classical Ferrari method.
    // Applying the substitution
    //     x = y - a/4
    // eliminates the cubic term and produces the depressed quartic
    //     y^4 + p*y^2 + q*y + r = 0.
    // Ferrari's method introduces an auxiliary variable m, obtained from
    // a resolvent cubic, and reduces the depressed quartic to two pairs
    // of quadratic solutions.

    // Shift used to eliminate the cubic term:
    //     x = y - shift.
    const double shift = a/4.0;
    const double shift_squared = shift*shift;

    // Coefficients of the depressed quartic
    //     y^4 + p*y^2 + q*y + r = 0.
    const double p = b - 6.0*shift_squared;
    const double q = c + (8.0*shift_squared - 2.0*b)*shift;
    const double r = d - c*shift + (b - 3.0*shift_squared)*shift_squared;

    // Auxiliary expression appearing in Ferrari's
    // completion-of-squares construction:
    //     h = p^2/4 - r.
    const double h = p*p/4.0 - r;

    // Solve the resolvent cubic
    //     m^3 + p*m^2 + h*m - q^2/8 = 0.
    // A nonnegative m is required because the final root formulas
    // contain sqrt(m/2).
    const double m = solve_cubic_one_real_root(p, h, -q*q/8.0);
    if (m < 0.0) return;

    // Ferrari's construction gives the coupling term
    //     coupling^2 = (m + p)*m + h.
    const double coupling_radicand = (m+p)*m + h;
    if (coupling_radicand < 0.0) return;
    const double sign_q =(q > 0.0) ? 1.0 : -1.0;
    const double coupling = sign_q*std::sqrt(coupling_radicand);

    // Common center of the two remaining radicands:
    const double radicand_center = -(m+p)/2.0;

    // Two pairs of possible depressed-quartic roots:
    //     y = sqrt(m/2) +/- sqrt(-(m+p)/2 - coupling).
    //     y = -sqrt(m/2) +/- sqrt(-(m+p)/2 + coupling).
    const double first_radicand = radicand_center - coupling;
    const double second_radicand = radicand_center + coupling;

    // Common square-root term in the four root expressions.
    const double sqrt_half_m = std::sqrt(m/2.0);

    if (first_radicand >= 0.0) {
        const double sqrt_first_radicand = std::sqrt(first_radicand);

        // Undo the translation x = y - shift.
        roots[n_roots++] = sqrt_half_m + sqrt_first_radicand - shift ;
        roots[n_roots++] = sqrt_half_m - sqrt_first_radicand - shift ;
    }

    if (second_radicand >= 0.0) {
        const double sqrt_second_radicand = std::sqrt(second_radicand);

        // Undo the translation x = y - shift.
        roots[n_roots++] = -sqrt_half_m + sqrt_second_radicand - shift ;
        roots[n_roots++] = -sqrt_half_m - sqrt_second_radicand - shift ;
    }

}

static void newton_raphson(T& d1, T& d2, T& d3, T s12, T s13, T s23, T m12, T m13, T m23) {
    // Three constraints:
    //     d1^2 + d2^2 - 2*d1*d2*m12 - s12 = 0.
    //     d1^2 + d3^2 - 2*d1*d3*m13 - s13 = 0.
    //     d2^2 + d3^2 - 2*d2*d3*m23 - s23 = 0.
    // Newton–Raphson linearizes these three nonlinear equations and solves
    //     J*delta = -residuals
    // for a depth correction delta at each iteration.

    for(int i=0;i<refinement_iterations;++i){
        const T d1_squared= d1*d1;
        const T d2_squared= d2*d2;
        const T d3_squared= d3*d3;
        const T r1 = d1_squared + d2_squared - 2.0*d1*d2*m12 - s12;
        const T r2 = d1_squared + d3_squared - 2.0*d1*d3*m13 - s13;
        const T r3 = d2_squared + d3_squared - 2.0*d2*d3*m23 - s23;
        if (std::abs(r1) + std::abs(r2) + std::abs(r3) < 1e-11) return;

        // After removing a common factor of 2, the Jacobian is
        //         [ j11  j12   0  ]
        //     2 * [ j21   0   j23 ]
        //         [  0   j32  j33 ],
        // where, for example,
        //     (1/2) * dr1/dd1 = d1 - d2*m12,
        //     (1/2) * dr1/dd2 = d2 - d1*m12.
        // The structural zeros allow the Newton step to be written
        // explicitly without constructing or inverting a matrix.
        const T j11 = d1 - d2*m12;
        const T j12 = d2 - d1*m12;
        const T j21 = d1 - d3*m13;
        const T j23 = d3 - d1*m13;
        const T j32 = d2 - d3*m23;
        const T j33 = d3 - d2*m23;

        // For the half-Jacobian above,
        //     det(J/2) = -(j11*j23*j32 + j12*j21*j33).
        // The factor 0.5 accounts for the common factor of 2 in the
        // full Jacobian. The remaining signs are incorporated directly
        // into the closed-form Newton updates below.
        const T inverse_determinant_factor = 0.5/(j11*j23*j32 + j12*j21*j33);

        // Apply the closed-form solution of J*delta = -residuals.
        // These expressions are the entries of the inverse Jacobian
        // multiplied by the residual vector.
        d1 += (-j23*j32*r1 - j12*j33*r2 + j12*j23*r3)*inverse_determinant_factor;
        d2 += (-j21*j33*r1 + j11*j33*r2 - j11*j23*r3)*inverse_determinant_factor;
        d3 +=(j21*j32*r1 - j11*j32*r2 - j12*j21*r3)*inverse_determinant_factor;
    }
}

public:
int solve(std::array<T,3> y1, std::array<T,3> y2, std::array<T,3> y3,
          std::array<T,3> x1, std::array<T,3> x2, std::array<T,3> x3,
          std::array<std::array<std::array<T,3>,3>,4>& Rs,
          std::array<std::array<T,3>,4>& Ts) const {

    // Normalize the rays
    const T y1_norm=std::sqrt(y1[0]*y1[0]+y1[1]*y1[1]+y1[2]*y1[2]);
    const T y2_norm=std::sqrt(y2[0]*y2[0]+y2[1]*y2[1]+y2[2]*y2[2]);
    const T y3_norm=std::sqrt(y3[0]*y3[0]+y3[1]*y3[1]+y3[2]*y3[2]);
    y1={y1[0]/y1_norm,y1[1]/y1_norm,y1[2]/y1_norm};
    y2={y2[0]/y2_norm,y2[1]/y2_norm,y2[2]/y2_norm};
    y3={y3[0]/y3_norm,y3[1]/y3_norm,y3[2]/y3_norm};

    // Get dot products
    T m12=y1[0]*y2[0]+y1[1]*y2[1]+y1[2]*y2[2];
    T m13=y1[0]*y3[0]+y1[1]*y3[1]+y1[2]*y3[2];
    T m23=y2[0]*y3[0]+y2[1]*y3[1]+y2[2]*y3[2];

    // Reindex if neccessary to keep m13 minimum:
    if (m12<=m23 && m12<=m13) { std::swap(x2,x3); std::swap(y2,y3); std::swap(m12,m13); } // Swap 2 and 3
    else if (m23<=m12 && m23<=m13) { std::swap(x1,x2); std::swap(y1,y2); std::swap(m13,m23); } // Swap 1 and 2

    // Reindex if neccessary to keep m23 maximum:
    if (m12>m23) { std::swap(x1,x3); std::swap(y1,y3); std::swap(m12,m23); } // Swap 1 and 3

    const std::array<T,3> d12{x1[0]-x2[0],x1[1]-x2[1],x1[2]-x2[2]}; // x1-x2
    const std::array<T,3> d13{x1[0]-x3[0],x1[1]-x3[1],x1[2]-x3[2]}; // x1-x3
    const std::array<T,3> d23{x2[0]-x3[0],x2[1]-x3[1],x2[2]-x3[2]}; // x2-x3

    // Get squared norms
    const T s12=d12[0]*d12[0]+d12[1]*d12[1]+d12[2]*d12[2];
    const T s13=d13[0]*d13[0]+d13[1]*d13[1]+d13[2]*d13[2];
    const T s23=d23[0]*d23[0]+d23[1]*d23[1]+d23[2]*d23[2];

    // Compute the coefficients
    const T s12_2=s12*s12, s13_2=s13*s13, s23_2=s23*s23;
    const T m12_2=m12*m12, m13_2=m13*m13, m23_2=m23*m23;
    const T s12s13=s12*s13, s12s23=s12*s23, s13s23=s13*s23;
    const T m12m23=m12*m23, m12m13m23=m12m23*m13;
    T coeff[5];
    coeff[0]=-s12_2+2*s12s13+2*s12s23-s13_2+4*s13s23*m12_2-2*s13s23-s23_2;
    coeff[1]=4*s12_2*m13-4*s12s13*m12m23-4*s12s13*m13-8*s12s23*m13+4*s13_2*m12m23-8*s13s23*m12_2*m13-4*s13s23*m12m23+4*s13s23*m13+4*s23_2*m13;
    coeff[2]=-4*s12_2*m13_2-2*s12_2+8*s12s13*m12m13m23+4*s12s13*m23_2+8*s12s23*m13_2+4*s12s23-4*s13_2*m12_2-4*s13_2*m23_2+2*s13_2+4*s13s23*m12_2+8*s13s23*m12m13m23-4*s23_2*m13_2-2*s23_2;
    coeff[3]=4*s12_2*m13-4*s12s13*m12m23-8*s12s13*m13*m23_2+4*s12s13*m13-8*s12s23*m13+4*s13_2*m12m23-4*s13s23*m12m23-4*s13s23*m13+4*s23_2*m13;
    coeff[4]=-s12_2+4*s12s13*m23_2-2*s12s13+2*s12s23-s13_2+2*s13s23-s23_2;

    const double c0_inv=1.0/coeff[0], a=coeff[1]*c0_inv, b=coeff[2]*c0_inv, c=coeff[3]*c0_inv, d=coeff[4]*c0_inv;
    double roots[4]{}; int n_roots=0;

    // Solve the quartic
    if (std::abs(a)>10.0) solve_quartic_FerrariLagrange(a,b,c,d,roots,n_roots);
    else solve_quartic_ClassicalFerrari(a,b,c,d,roots,n_roots);

    // Compute matrix X^(-1)
    const std::array<T,3> n{d12[1]*d13[2]-d12[2]*d13[1],
                            d12[2]*d13[0]-d12[0]*d13[2],
                            d12[0]*d13[1]-d12[1]*d13[0]};
    const double one_over_det=1.0/(d12[0]*(d13[1]*n[2]-n[1]*d13[2])-d13[0]*(d12[1]*n[2]-n[1]*d12[2])+n[0]*(d12[1]*d13[2]-d13[1]*d12[2]));
    const std::array<std::array<T,3>,3> X_inv{{
                                             {{(d13[1]*n[2]-n[1]*d13[2])*one_over_det,
                                               (n[0]*d13[2]-d13[0]*n[2])*one_over_det,
                                               (d13[0]*n[1]-n[0]*d13[1])*one_over_det}},
                                             {{(n[1]*d12[2]-d12[1]*n[2])*one_over_det,
                                               (d12[0]*n[2]-n[0]*d12[2])*one_over_det,
                                               (n[0]*d12[1]-d12[0]*n[1])*one_over_det}},
                                             {{(d12[1]*d13[2]-d13[1]*d12[2])*one_over_det,
                                               (d13[0]*d12[2]-d12[0]*d13[2])*one_over_det,
                                               (d12[0]*d13[1]-d13[0]*d12[1])*one_over_det}}
    }};

    // Recover R and t
    int valid=0;
    for (int i=0; i<n_roots && valid<4; ++i) {
        const T x=roots[i];
        if (x<=0) continue;
        const T den=2*s13*(m12*x-m23);
        if (den==0) continue;
        const T num=(s13+s23-s12)*x*x+2*m13*(s12-s23)*x+(s23-s12-s13);
        if ((num<0&&den>0)||(num>0&&den<0)) continue;
        const T y=num/den;
        const T d3_den=y*y-2*y*m23+1;
        if (d3_den<=0) continue;
        T d3=std::sqrt(s23/d3_den), d1=x*d3, d2=y*d3;

        // Refine the depths
        newton_raphson(d1,d2,d3,s12,s13,s23,m12,m13,m23);

        // Compute matrix Y
        const std::array<T,3> yd12{y1[0]*d1-y2[0]*d2,y1[1]*d1-y2[1]*d2,y1[2]*d1-y2[2]*d2};
        const std::array<T,3> yd13{y1[0]*d1-y3[0]*d3,y1[1]*d1-y3[1]*d3,y1[2]*d1-y3[2]*d3};
        const std::array<T,3> yn{yd12[1]*yd13[2]-yd12[2]*yd13[1],
                                 yd12[2]*yd13[0]-yd12[0]*yd13[2],
                                 yd12[0]*yd13[1]-yd12[1]*yd13[0]};
        const std::array<std::array<T,3>,3> Y{{
            {{yd12[0],yd13[0],yn[0]}},
            {{yd12[1],yd13[1],yn[1]}},
            {{yd12[2],yd13[2],yn[2]}}
        }};

        // Rs[valid]=Y*X_inv
        for (std::size_t row=0;row<3;++row)
            for (std::size_t col=0;col<3;++col)
                Rs[valid][row][col]=Y[row][0]*X_inv[0][col]
                                    +Y[row][1]*X_inv[1][col]
                                    +Y[row][2]*X_inv[2][col];

        const std::array<T,3> rotated_x1{
            Rs[valid][0][0]*x1[0]+Rs[valid][0][1]*x1[1]+Rs[valid][0][2]*x1[2],
            Rs[valid][1][0]*x1[0]+Rs[valid][1][1]*x1[1]+Rs[valid][1][2]*x1[2],
            Rs[valid][2][0]*x1[0]+Rs[valid][2][1]*x1[1]+Rs[valid][2][2]*x1[2]};

        // Ts[valid]=(y1*d1 - Rs[valid]*x1)
        Ts[valid]={y1[0]*d1-rotated_x1[0],y1[1]*d1-rotated_x1[1],y1[2]*d1-rotated_x1[2]};
        ++valid;
    }
    return valid;
}
};

} // namespace p3p
