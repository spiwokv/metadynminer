#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix hills1(std::vector<double> cv1, std::vector<double> cv2, double width1, double width2, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i, cv2i;
  int ni, nj;
  double z;
  double g[2*n][2*n];
  NumericMatrix v(n, n);
  for (int i = 0; i < 2*n; i++) {
    for (int j = 0; j < 2*n; j++) {
      g[i][j] = 0.0;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      z = exp(-double(i)*double(i)/2.0/width1/width1-double(j)*double(j)/2.0/width2/width2);
      g[i][j] = z;
      if(j>0) g[i][2*n-j] = z;
      if(i>0) g[2*n-i][j] = z;
      if((i>0)&&(j>0)) g[2*n-i][2*n-j] = z;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j) = 0.0;
    }
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    cv2i = int(cv2[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=2*n;
      for (int j = 0; j < n; j++) {
        nj=j-cv2i;
        if(j<cv2i) nj+=2*n;
        v(i, j) -= heights[icv]*g[ni][nj];
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills1p1(std::vector<double> cv1, std::vector<double> cv2, double width1, double width2, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i, cv2i;
  int ni, nj;
  double z;
  double g[n][2*n];
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2*n; j++) {
      g[i][j] = 0.0;
    }
  }
  for (int i = 0; i < n/2; i++) {
    for (int j = 0; j < n; j++) {
      z = exp(-double(i)*double(i)/2.0/width1/width1-double(j)*double(j)/2.0/width2/width2);
      g[i][j] = z;
      if(j>0) g[i][2*n-j] = z;
      if(i>0) g[n-i][j] = z;
      if((i>0)&&(j>0)) g[n-i][2*n-j] = z;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j) = 0.0;
    }
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    cv2i = int(cv2[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=n;
      for (int j = 0; j < n; j++) {
        nj=j-cv2i;
        if(j<cv2i) nj+=2*n;
        v(i, j) -= heights[icv]*g[ni][nj];
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills1p2(std::vector<double> cv1, std::vector<double> cv2, double width1, double width2, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i, cv2i;
  int ni, nj;
  double z;
  double g[2*n][n];
  NumericMatrix v(n, n);
  for (int i = 0; i < 2*n; i++) {
    for (int j = 0; j < n; j++) {
      g[i][j] = 0.0;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n/2; j++) {
      z = exp(-double(i)*double(i)/2.0/width1/width1-double(j)*double(j)/2.0/width2/width2);
      g[i][j] = z;
      if(j>0) g[i][n-j] = z;
      if(i>0) g[2*n-i][j] = z;
      if((i>0)&&(j>0)) g[2*n-i][n-j] = z;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j) = 0.0;
    }
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    cv2i = int(cv2[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=2*n;
      for (int j = 0; j < n; j++) {
        nj=j-cv2i;
        if(j<cv2i) nj+=n;
        v(i, j) -= heights[icv]*g[ni][nj];
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills1p12(std::vector<double> cv1, std::vector<double> cv2, double width1, double width2, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i, cv2i;
  int ni, nj;
  double z;
  double g[n][n];
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      g[i][j] = 0.0;
    }
  }
  for (int i = 0; i < n/2; i++) {
    for (int j = 0; j < n/2; j++) {
      z = exp(-double(i)*double(i)/2.0/width1/width1-double(j)*double(j)/2.0/width2/width2);
      g[i][j] = z;
      if(j>0) g[i][n-j] = z;
      if(i>0) g[n-i][j] = z;
      if((i>0)&&(j>0)) g[n-i][n-j] = z;
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j) = 0.0;
    }
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    cv2i = int(cv2[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=n;
      for (int j = 0; j < n; j++) {
        nj=j-cv2i;
        if(j<cv2i) nj+=n;
        v(i, j) -= heights[icv]*g[ni][nj];
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills2(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j)=0.0;
    }
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        dcv1 = cv1[icv]-double(i);
        dcv2 = cv2[icv]-double(j);
        v(i, j) -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]-dcv2*dcv2/2.0/width2[icv]/width2[icv]);
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills2p1(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j)=0.0;
    }
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        dcv1 = cv1[icv]-double(i);
        dcv2 = cv2[icv]-double(j);
        if(dcv1 >  double(n)/2.0) dcv1 -= double(n);
        if(dcv1 < -double(n)/2.0) dcv1 += double(n);
        v(i, j) -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]-dcv2*dcv2/2.0/width2[icv]/width2[icv]);
      }
    }
  }
  return v;
}


// [[Rcpp::export]]
NumericMatrix hills2p2(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j)=0.0;
    }
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        dcv1 = cv1[icv]-double(i);
        dcv2 = cv2[icv]-double(j);
        if(dcv2 >  double(n)/2.0) dcv2 -= double(n);
        if(dcv2 < -double(n)/2.0) dcv2 += double(n);
        v(i, j) -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]-dcv2*dcv2/2.0/width2[icv]/width2[icv]);
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
NumericMatrix hills2p12(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  NumericMatrix v(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      v(i, j)=0.0;
    }
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        dcv1 = cv1[icv]-double(i);
        dcv2 = cv2[icv]-double(j);
        if(dcv1 >  double(n)/2.0) dcv1 -= double(n);
        if(dcv1 < -double(n)/2.0) dcv1 += double(n);
        if(dcv2 >  double(n)/2.0) dcv2 -= double(n);
        if(dcv2 < -double(n)/2.0) dcv2 += double(n);
        v(i, j) -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]-dcv2*dcv2/2.0/width2[icv]/width2[icv]);
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
std::vector<double> hills1d1(std::vector<double> cv1, double width1, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i;
  int ni;
  double z;
  double g[2*n];
  std::vector<double> v(n);
  for (int i = 0; i < 2*n; i++) {
    g[i]=0.0;
  }
  for (int i = 0; i < n; i++) {
    z = exp(-double(i)*double(i)/2.0/width1/width1);
    g[i] = z;
    if(i>0) g[2*n-i] = z;
  }
  for (int i = 0; i < n; i++) {
    v[i] = 0.0;
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=2*n;
      v[i] -= heights[icv]*g[ni];
    }
  }
  return v;
}

// [[Rcpp::export]]
std::vector<double> hills1d1p(std::vector<double> cv1, double width1, std::vector<double> heights, int n, int tmin, int tmax) {
  int cv1i;
  int ni;
  double z;
  double g[n];
  std::vector<double> v(n);
  for (int i = 0; i < n; i++) {
    g[i]=0.0;
  }
  for (int i = 0; i < n/2; i++) {
    z = exp(-double(i)*double(i)/2.0/width1/width1);
    g[i] = z;
    if(i>0) g[n-i] = z;
  }
  for (int i = 0; i < n; i++) {
    v[i] = 0.0;
  }
  for (int icv=tmin; icv<=tmax; icv++) {
    cv1i = int(cv1[icv]+0.5);
    for (int i = 0; i < n; i++) {
      ni=i-cv1i;
      if(i<cv1i) ni+=n;
      v[i] -= heights[icv]*g[ni];
    }
  }
  return v;
}

// [[Rcpp::export]]
std::vector<double> hills1d2(std::vector<double> cv1, std::vector<double> width1, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  std::vector<double> v(n);
  for (int i = 0; i < n; i++) {
    v[i]=0.0;
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      dcv1 = cv1[icv]-double(i);
      v[i] -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]);
    }
  }
  return v;
}

// [[Rcpp::export]]
std::vector<double> hills1d2p(std::vector<double> cv1, std::vector<double> width1, std::vector<double> heights, int n, int tmin, int tmax) {
  double dcv1;
  std::vector<double> v(n);
  for (int i = 0; i < n; i++) {
    v[i]=0.0;
  }
  for (int icv=tmin; icv <= tmax; icv++) {
    for (int i = 0; i < n; i++) {
      dcv1 = cv1[icv]-double(i);
      if(dcv1 >  double(n)/2.0) dcv1 -= double(n);
      if(dcv1 < -double(n)/2.0) dcv1 += double(n);
      v[i] -= heights[icv]*exp(-dcv1*dcv1/2.0/width1[icv]/width1[icv]);
    }
  }
  return v;
}


// [[Rcpp::export]]
std::vector<double> fe2d(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, double x, double y, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    dcv2 = cv2[i]-y;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]-dcv2*dcv2/2.0/width2[i]/width2[i]);
    vo[i]=v;
  }
  return vo;
}

// [[Rcpp::export]]
std::vector<double> fe2dp1(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, double x, double y, double p1, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    if(dcv1 >  p1/2.0) dcv1 -= p1;
    if(dcv1 < -p1/2.0) dcv1 += p1;
    dcv2 = cv2[i]-y;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]-dcv2*dcv2/2.0/width2[i]/width2[i]);
    vo[i]=v;
  }
  return vo;
}

// [[Rcpp::export]]
std::vector<double> fe2dp2(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, double x, double y, double p2, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    dcv2 = cv2[i]-y;
    if(dcv2 >  p2/2.0) dcv2 -= p2;
    if(dcv2 < -p2/2.0) dcv2 += p2;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]-dcv2*dcv2/2.0/width2[i]/width2[i]);
    vo[i]=v;
  }
  return vo;
}

// [[Rcpp::export]]
std::vector<double> fe2dp12(std::vector<double> cv1, std::vector<double> cv2, std::vector<double> width1, std::vector<double> width2, std::vector<double> heights, double x, double y, double p1, double p2, int tmin, int tmax) {
  double dcv1;
  double dcv2;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    if(dcv1 >  p1/2.0) dcv1 -= p1;
    if(dcv1 < -p1/2.0) dcv1 += p1;
    dcv2 = cv2[i]-y;
    if(dcv2 >  p2/2.0) dcv2 -= p2;
    if(dcv2 < -p2/2.0) dcv2 += p2;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]-dcv2*dcv2/2.0/width2[i]/width2[i]);
    vo[i]=v;
  }
  return vo;
}

// [[Rcpp::export]]
std::vector<double> fe1d(std::vector<double> cv1, std::vector<double> width1, std::vector<double> heights, double x, int tmin, int tmax) {
  double dcv1;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]);
    vo[i]=v;
  }
  return vo;
}

// [[Rcpp::export]]
std::vector<double> fe1dp(std::vector<double> cv1, std::vector<double> width1, std::vector<double> heights, double x, double p1, int tmin, int tmax) {
  double dcv1;
  double v;
  std::vector<double> vo(tmax-tmin+1);
  v = 0.0;
  for (int i=tmin; i <= tmax; i++) {
    dcv1 = cv1[i]-x;
    if(dcv1 >  p1/2.0) dcv1 -= p1;
    if(dcv1 < -p1/2.0) dcv1 += p1;
    v -= heights[i]*exp(-dcv1*dcv1/2.0/width1[i]/width1[i]);
    vo[i]=v;
  }
  return vo;
}
