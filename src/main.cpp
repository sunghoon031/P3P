
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include "p3p.h"


std::array<std::array<double,3>,3> axis_angle(const std::array<double,3>& axis, double angle) {
    const double axis_norm=std::sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    const std::array<double,3> a{axis[0]/axis_norm,axis[1]/axis_norm,axis[2]/axis_norm};
    const double c=std::cos(angle),s=std::sin(angle),v=1-c;
    return {{{{c+a[0]*a[0]*v,a[0]*a[1]*v-a[2]*s,a[0]*a[2]*v+a[1]*s}},
             {{a[1]*a[0]*v+a[2]*s,c+a[1]*a[1]*v,a[1]*a[2]*v-a[0]*s}},
             {{a[2]*a[0]*v-a[1]*s,a[2]*a[1]*v+a[0]*s,c+a[2]*a[2]*v}}}};
}

double pose_error_l1(const std::array<std::array<double,3>,3>& R_est,
                     const std::array<double,3>& t_est,
                     const std::array<std::array<double,3>,3>& R_gt,
                     const std::array<double,3>& t_gt) {
    double error = 0.0;
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t col = 0; col < 3; ++col)
            error += std::abs(R_gt[row][col] - R_est[row][col]);
    for (std::size_t i = 0; i < 3; ++i)
        error += std::abs(t_gt[i] - t_est[i]);
    return error;
}

double determinant(const std::array<std::array<double,3>,3>& R) {
    return R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
    - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
        + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
}

double orthogonality_error_l1(const std::array<std::array<double,3>,3>& R) {
    double error=0.0;
    for (std::size_t row=0;row<3;++row)
        for (std::size_t col=0;col<3;++col) {
            const double identity=(row==col)?1.0:0.0;
            const double value=R[0][row]*R[0][col]
                                 +R[1][row]*R[1][col]
                                 +R[2][row]*R[2][col];
            error+=std::abs(value-identity);
        }
    return error;
}

double quaternion_norm_from_rotation_matrix(const std::array<std::array<double,3>,3>& R) {
    // Standard trace/diagonal conversion, equivalent to constructing a
    // quaternion directly from the matrix and then taking its Euclidean norm.
    double w=0.0,x=0.0,y=0.0,z=0.0;
    const double trace=R[0][0]+R[1][1]+R[2][2];
    if (trace>0.0) {
        const double s=2.0*std::sqrt(std::max(0.0,trace+1.0));
        if (s==0.0) return 0.0;
        w=0.25*s;
        x=(R[2][1]-R[1][2])/s;
        y=(R[0][2]-R[2][0])/s;
        z=(R[1][0]-R[0][1])/s;
    } else if (R[0][0]>R[1][1] && R[0][0]>R[2][2]) {
        const double s=2.0*std::sqrt(std::max(0.0,1.0+R[0][0]-R[1][1]-R[2][2]));
        if (s==0.0) return 0.0;
        w=(R[2][1]-R[1][2])/s;
        x=0.25*s;
        y=(R[0][1]+R[1][0])/s;
        z=(R[0][2]+R[2][0])/s;
    } else if (R[1][1]>R[2][2]) {
        const double s=2.0*std::sqrt(std::max(0.0,1.0+R[1][1]-R[0][0]-R[2][2]));
        if (s==0.0) return 0.0;
        w=(R[0][2]-R[2][0])/s;
        x=(R[0][1]+R[1][0])/s;
        y=0.25*s;
        z=(R[1][2]+R[2][1])/s;
    } else {
        const double s=2.0*std::sqrt(std::max(0.0,1.0+R[2][2]-R[0][0]-R[1][1]));
        if (s==0.0) return 0.0;
        w=(R[1][0]-R[0][1])/s;
        x=(R[0][2]+R[2][0])/s;
        y=(R[1][2]+R[2][1])/s;
        z=0.25*s;
    }
    return std::sqrt(w*w+x*x+y*y+z*z);
}

struct RotationValidity {
    double determinant_error{};
    double orthogonality_error{};
    double quaternion_norm_error{};
    bool valid{};
};

RotationValidity check_rotation_validity(const std::array<std::array<double,3>,3>& R) {
    RotationValidity result;
    result.determinant_error=std::abs(determinant(R)-1.0);
    result.orthogonality_error=orthogonality_error_l1(R);
    result.quaternion_norm_error=std::abs(1.0-quaternion_norm_from_rotation_matrix(R));
    result.valid=result.determinant_error<1e-6
                   && result.orthogonality_error<1e-6
                   && result.quaternion_norm_error<1e-5;
    return result;
}


int main() {
    std::random_device random_device;
    std::uniform_int_distribution<unsigned int> seed_distribution(0,99);
    const unsigned int seed=seed_distribution(random_device);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> image(-1.0,1.0),depth(0.1,10.0),unit(-1.0,1.0),angle(-3.14159,3.14159);
    std::array<std::array<double,3>,3> rays{},camera_points{},world_points{};
    for (std::size_t i=0;i<3;++i) {
        const double u=image(rng),v=image(rng);
        const double ray_norm=std::sqrt(u*u+v*v+1.0);
        rays[i]={u/ray_norm,v/ray_norm,1.0/ray_norm};
        const double point_depth=depth(rng);
        camera_points[i]={point_depth*rays[i][0],point_depth*rays[i][1],point_depth*rays[i][2]};
    }
    const std::array<std::array<double,3>,3> R_gt=axis_angle({unit(rng),unit(rng),unit(rng)},angle(rng));
    const std::array<double,3> raw_translation{unit(rng),unit(rng),unit(rng)};
    const double translation_norm=std::sqrt(raw_translation[0]*raw_translation[0]
                                              +raw_translation[1]*raw_translation[1]
                                              +raw_translation[2]*raw_translation[2]);
    const std::array<double,3> t_gt{raw_translation[0]/translation_norm,
                                     raw_translation[1]/translation_norm,
                                     raw_translation[2]/translation_norm};
    const std::array<std::array<double,3>,3> R_gt_T{{
        {{R_gt[0][0],R_gt[1][0],R_gt[2][0]}},
        {{R_gt[0][1],R_gt[1][1],R_gt[2][1]}},
        {{R_gt[0][2],R_gt[1][2],R_gt[2][2]}}
    }};
    for (std::size_t i=0;i<3;++i) {
        const std::array<double,3> shifted{camera_points[i][0]-t_gt[0],
                                            camera_points[i][1]-t_gt[1],
                                            camera_points[i][2]-t_gt[2]};
        world_points[i]={R_gt_T[0][0]*shifted[0]+R_gt_T[0][1]*shifted[1]+R_gt_T[0][2]*shifted[2],
                           R_gt_T[1][0]*shifted[0]+R_gt_T[1][1]*shifted[1]+R_gt_T[1][2]*shifted[2],
                           R_gt_T[2][0]*shifted[0]+R_gt_T[2][1]*shifted[1]+R_gt_T[2][2]*shifted[2]};
    }

    std::array<std::array<std::array<double,3>,3>,4> Rs{};
    std::array<std::array<double,3>,4> ts{};
    const p3p::P3P<double,5> solver;
    const int count=solver.solve(rays[0],rays[1],rays[2],
                                   world_points[0],world_points[1],world_points[2],
                                   Rs,ts);
    std::cout<<std::fixed<<std::setprecision(10)<<"Synthetic P3P instance (seed "<<seed<<")\nNumber of Solutions: "<<count<<"\n";
    bool found_ground_truth_solution=false;
    for (int k=0;k<count;++k) {
        const double error=pose_error_l1(Rs[k],ts[k],R_gt,t_gt);
        const RotationValidity validity=check_rotation_validity(Rs[k]);
        std::cout<<"Solution "<<k<<":\n"
                  <<"  pose error (L1)       = "<<error<<'\n'
                  <<"  |det(R) - 1|          = "<<validity.determinant_error<<'\n'
                  <<"  ||R^T R - I||_L1      = "<<validity.orthogonality_error<<'\n'
                  <<"  |1 - quaternion norm| = "<<validity.quaternion_norm_error<<'\n'
                  <<"  valid rotation        = "<<(validity.valid?"yes":"no")<<'\n';
        if (validity.valid && error<1e-6) found_ground_truth_solution=true;
    }
    // The paper's ground-truth classification requires a valid returned pose
    // whose combined L1 pose error is below 1e-6.
    const bool passed=found_ground_truth_solution;
    std::cout<<"GT Verification: "<<(passed?"PASS":"FAIL")<<'\n';
    return passed?0:1;
}
