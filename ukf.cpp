
#include "Eigen/Dense"
#include "ukf.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF(){
  use_laser_=true;
  use_radar_=true;
  n_x_=5;
  x_=VectorXd(n_x_);
  P_=MatrixXd(n_x_,n_x_);
  P_<< 1,0,0,0,0,
      0,1,0,0,0,
      0,0,1,0,0,
      0,0,0,0.0225,0,
      0,0,0,0,0.0225;
  std_yawdd_=0.5;
  std_a_=0.8;
  std_laspx_=0.15;
  std_laspy_=0.15;
  std_radr_=0.3;
  std_radphi_=0.03;
  std_radrd_=0.3;

  n_aug_=n_x_+2;
  lambda_=3-n_aug_;
  time_us_=0.0;
  weights_=VectorXd(2*n_aug_ +1);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0)=lambda_/(lambda_+n_aug_);

  Xsig_pred_=MatrixXd(n_x_,2*n_aug_+1);
  R_radar_=MatrixXd(3,3);
  R_radar_<<std_radr_*std_radr_,0,0,
          0,std_radphi_*std_radphi_,0,
          0,0,std_radrd_*std_radrd_;

  R_lidar_= MatrixXd(2,2);
  R_lidar_<<std_laspx_*std_laspx_,0,
            0,std_laspy_*std_laspy_;

  NIS_laser_=0;
  NIS_radar_=0;
}

UKF::~UKF(){}
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_==MeasurementPackage::RADAR)
    {
      double rho= meas_package.raw_measurements_[0];
      double phi=meas_package.raw_measurements_[1];
      double rho_dot=meas_package.raw_measurements_[2];
      double x=rho*cos(phi);
      double y=rho*sin(phi);
      double vx=rho_dot*cos(phi);
      double vy=rho_dot*sin(phi);
      double v=sqrt(vx*vx+vy*vy);
      x_<< x,
           y,
           v,
           0,
           0;
    }
    else
    {
      x_<<meas_package.raw_measurements_[0],
         meas_package.raw_measurements_[1],
         0,
         0,
         0;
    }
    time_us_=meas_package.timestamp_;
    is_initialized_=true;
    return;    
  }
  double dt=(meas_package.timestamp_-time_us_)/1000000.0;
  time_us_=meas_package.timestamp_;
  Prediction(dt);

  if(meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  if(meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }

}
void UKF::Prediction(double delta_t)
{
  VectorXd x_aug=VectorXd(n_aug_);
  MatrixXd P_aug=MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);
  x_aug.head(5)=x_;
  x_aug(5)=0;
  x_aug(6)=0;

  P_aug.fill(0);
  P_aug.topLeftCorner(5,5)=P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6)=std_yawdd_*std_yawdd_;

  MatrixXd L=P_aug.llt().matrixL();
  Xsig_aug.col(0)=x_aug;
  for(int i=0;i<n_aug_;i++)
  {
    Xsig_aug.col(i+1)=x_aug+sqrt(lambda_+n_aug_)*L.col(i);
    Xsig_aug.col(i+1+n_aug_)=x_aug-sqrt(lambda_+n_aug_)*L.col(i);
  }
  for(int i=0;i<2*n_aug_+1;i++)
  {
    double p_x=Xsig_aug(0,i);
    double p_y=Xsig_aug(1,i);
    double v=Xsig_aug(2,i);
    double yaw=Xsig_aug(3,i);
    double yawd= Xsig_aug(4,i);
    double nu_a=Xsig_aug(5,i);
    double nu_yawd=Xsig_aug(6,i);
    
    double px_p,py_p,v_p,yaw_p,yawd_p;
    if(fabs(yawd)>0.001)
    {
      px_p=p_x+v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw));
      py_p=p_y+v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw));
    }
    else
    {
      px_p=p_x+v*delta_t*cos(yaw);
      py_p=p_y+v*delta_t*sin(yaw);
    }
    v_p=v;
    yaw_p=yaw+yawd*delta_t;
    yawd_p=yawd;

    px_p=px_p+0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p=py_p+0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p=v_p+nu_a*delta_t;
    yaw_p=yaw_p+0.5*nu_yawd*delta_t*delta_t;
    yawd_p=yawd_p+nu_yawd*delta_t;

    Xsig_pred_(0,i)=px_p;
    Xsig_pred_(1,i)=py_p;
    Xsig_pred_(2,i)=v_p;
    Xsig_pred_(3,i)=yaw_p;
    Xsig_pred_(4,i)=yawd_p;
  }
  x_.fill(0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    x_=x_+weights_(i)*Xsig_pred_.col(i);
  }
  P_.fill(0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2.*-M_PI;
    P_=P_+weights_(i)*x_diff*x_diff.transpose();
  }

}
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  VectorXd z_=meas_package.raw_measurements_;
  int n_z_=2;
  MatrixXd Zsig=MatrixXd(n_z_,2*n_aug_+1);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    Zsig(0,i)=Xsig_pred_(0,i);
    Zsig(1,i)=Xsig_pred_(1,i);
  }
  VectorXd z_pred_=VectorXd(n_z_);
  z_pred_.fill(0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    z_pred_=z_pred_+weights_(i)*Zsig.col(i);
  }
  MatrixXd S=MatrixXd(n_z_,n_z_);
  S.fill(0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    VectorXd z_diff=Zsig.col(i)-z_pred_;
    S=S+weights_(i)*z_diff*z_diff.transpose();
  }

  S=S+R_lidar_;

  MatrixXd Tc=MatrixXd(n_x_,n_z_);
  Tc.fill(0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    VectorXd z_diff=Zsig.col(i)-z_pred_;
    Tc=Tc+weights_(i)*x_diff*z_diff.transpose();
  }
  MatrixXd K=Tc*S.inverse();
  VectorXd z_diff=z_-z_pred_;
  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();
  NIS_laser_=z_diff.transpose()*S.inverse()*z_diff;
}
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  VectorXd z_=meas_package.raw_measurements_;
  int n_z_=3;
  MatrixXd Zsig=MatrixXd(n_z_,2*n_aug_+1);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    double p_x=Xsig_pred_(0,i);
    double p_y=Xsig_pred_(1,i);
    double v=Xsig_pred_(2,i);
    double yaw=Xsig_pred_(3,i);
    double yawd=Xsig_pred_(4,i);
    double vx=v*cos(yaw);
    double vy=v*sin(yaw);

    Zsig(0,i)=sqrt(p_x*p_x+p_y*p_y);
    Zsig(1,i)=atan2(p_y,p_x);
    Zsig(2,i)=(p_x*vx+p_y*vy)/(sqrt(p_x*p_x+p_y*p_y));
  }
  VectorXd z_pred_=VectorXd(n_z_);
  z_pred_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    z_pred_=z_pred_+weights_(i)*Zsig.col(i);
  }

  MatrixXd S=MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    VectorXd z_diff=Zsig.col(i)-z_pred_;
    while(z_diff(1)>M_PI) 
    {
      z_diff(1)-=2.*M_PI;
    }  
    while(z_diff(1)<-M_PI) 
    {
      z_diff(1)+=2.*M_PI;
    }  
    S=S+weights_(i)*z_diff*z_diff.transpose();
  }
  S=S+R_radar_;

  MatrixXd Tc=MatrixXd(n_x_,n_z_);
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    VectorXd x_diff=Xsig_pred_.col(i)-x_;
    while(x_diff(3)>M_PI) 
    {
      x_diff(3)-=2.*M_PI;
    }  
    while(x_diff(3)<-M_PI) 
    {
      x_diff(3)+=2*M_PI;
    }

    VectorXd z_diff=Zsig.col(i)-z_pred_;
    while(z_diff(1)>M_PI) 
    {
      z_diff(1)-=2.*M_PI;
    }
    while(z_diff(1)<-M_PI) 
    {
      z_diff(1)+=2.*M_PI;
    }
    Tc=Tc+weights_(i)*x_diff*z_diff.transpose();
  }
  MatrixXd K=Tc*S.inverse();
  VectorXd z_diff=z_-z_pred_;
  while(z_diff(1)>M_PI) 
  {
    z_diff(1)-=2.*M_PI;
  }  
  while(z_diff(1)<-M_PI) 
  {
    z_diff(1)+=2.*M_PI;
  }
  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();
  NIS_radar_=z_diff.transpose()*S.inverse()*z_diff;
}






