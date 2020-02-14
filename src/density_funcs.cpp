// Functions to evaluate the DDM PDF for specific criteria

#include "funcs.h"




//////////                                                            //////////
///////////////////////////////// Small Time ///////////////////////////////////
//////////                                                            //////////

// Foster style with 2017 sum
double fs_Fos_17(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_eps_17(t, a, w, eps / exp(mult)));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_eps_17(t, a, w, eps / exp(mult)));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_eps_17(t, a, w, eps / mult);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_eps_17(t, a, w, eps / mult);
    }
  }
}


// Foster style with 2014 sum
double fs_Fos_14(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_eps_14(t, a, w, eps / exp(mult)));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_eps_14(t, a, w, eps / exp(mult)));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_eps_14(t, a, w, eps / mult);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_eps_14(t, a, w, eps / mult);
    }
  }
}


// Kesselmeier and co style with 2017 sum
double fs_Kes_17(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Kes(t / (a * a), w, eps);
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_2017(t, a, w, ks));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_2017(t, a, w, ks));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_2017(t, a, w, ks);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_2017(t, a, w, ks);
    }
  }
}


// Kesselmeier and co style with 2014 sum
double fs_Kes_14(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Kes(t / (a * a), w, eps);
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_2014(t, a, w, ks));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_2014(t, a, w, ks));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_2014(t, a, w, ks);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_2014(t, a, w, ks);
    }
  }
}


// Navarro style with 2017 sum
double fs_Nav_17(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Nav(t / (a * a), eps);
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_2017(t, a, w, ks));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_2017(t, a, w, ks));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_2017(t, a, w, ks);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_2017(t, a, w, ks);
    }
  }
}


// Navarro style with 2014 sum
double fs_Nav_14(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Nav(t / (a * a), eps);
  double mult;
  if (log_prob) { // log
    double log_a = log(a);
    double log_t3 = 1.5 * log(t);
    if (sv < SV_THRESH) { // no sv
      mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
      return mult + log(small_sum_2014(t, a, w, ks));
    } else { // sv
      double log_svt = 0.5 * log(1 + sv * t);
      mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
            - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
      return mult + log(small_sum_2014(t, a, w, ks));
    }
  } else { // no log
    if (sv < SV_THRESH) { // no sv
      mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
      return mult * small_sum_2014(t, a, w, ks);
    } else { // sv
      mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
            / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
      return mult * small_sum_2014(t, a, w, ks);
    }
  }
}



//////////                                                            //////////
///////////////////////////////// Large Time ///////////////////////////////////
//////////                                                            //////////

// Navarro style with large sum
double fl_Nav_09(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  // note: sv is not used because there is only a constant drift rate variant
  int kl = kl_Nav(t / (a * a), eps);
  double mult;
  if (log_prob) { // log
    double log_a2 = 2 * log(a);
    mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
    return mult + log(large_sum_Nav(t, a, w, kl));
  } else { // no log
    mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
    return mult * large_sum_Nav(t, a, w, kl);
  }
}



//////////                                                            //////////
///////////////////////////////// Both Times ///////////////////////////////////
//////////                                                            //////////

// ks = Kesselmeier and co, 2017 style sum approximation, kl = Navarro
double fb_Kes_17(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Kes(t / (a * a), w, eps);
  int kl = kl_Nav(t / (a * a), eps);
  double mult;
  if (ks < kl) { // small time appx is more efficient
    if (log_prob) { // log
      double log_a = log(a);
      double log_t3 = 1.5 * log(t);
      if (sv < SV_THRESH) { // no sv
        mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
        return mult + log(small_sum_2017(t, a, w, ks));
      } else { // sv
        double log_svt = 0.5 * log(1 + sv * t);
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
              - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        return mult + log(small_sum_2017(t, a, w, ks));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        return mult * small_sum_2017(t, a, w, ks);
      } else { // sv
        mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
              / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
        return mult * small_sum_2017(t, a, w, ks);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      double log_a2 = 2 * log(a);
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl));
    } else { // no log
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl);
    }
  }
}


// ks = Kesselmeier and co, 2014 style sum approximation, kl = Navarro
double fb_Kes_14(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Kes(t / (a * a), w, eps);
  int kl = kl_Nav(t / (a * a), eps);
  double mult;
  if (ks < kl) { // small time appx is more efficient
    if (log_prob) { // log
      double log_a = log(a);
      double log_t3 = 1.5 * log(t);
      if (sv < SV_THRESH) { // no sv
        mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
        return mult + log(small_sum_2014(t, a, w, ks));
      } else { // sv
        double log_svt = 0.5 * log(1 + sv * t);
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
              - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        return mult + log(small_sum_2014(t, a, w, ks));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        return mult * small_sum_2014(t, a, w, ks);
      } else { // sv
        mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
              / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
        return mult * small_sum_2014(t, a, w, ks);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      double log_a2 = 2 * log(a);
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl));
    } else { // no log
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl);
    }
  }
}


// ks = Navarro, 2017 style sum approximation, kl = Navarro
double fb_Nav_17(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Nav(t / (a * a), eps);
  int kl = kl_Nav(t / (a * a), eps);
  double mult;
  if (ks < kl) { // small time appx is more efficient
    if (log_prob) { // log
      double log_a = log(a);
      double log_t3 = 1.5 * log(t);
      if (sv < SV_THRESH) { // no sv
        mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
        return mult + log(small_sum_2017(t, a, w, ks));
      } else { // sv
        double log_svt = 0.5 * log(1 + sv * t);
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
              - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        return mult + log(small_sum_2017(t, a, w, ks));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        return mult * small_sum_2017(t, a, w, ks);
      } else { // sv
        mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
              / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
        return mult * small_sum_2017(t, a, w, ks);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      double log_a2 = 2 * log(a);
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl));
    } else { // no log
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl);
    }
  }
}


// ks = Navarro, 2014 style sum approximation, kl = Navarro
double fb_Nav_14(const double& t, const double& a, const double& v,
                 const double& w, const double& sv, const bool& log_prob,
                 const double& eps)
{
  int ks = ks_Nav(t / (a * a), eps);
  int kl = kl_Nav(t / (a * a), eps);
  double mult;
  if (ks < kl) { // small time appx is more efficient
    if (log_prob) { // log
      double log_a = log(a);
      double log_t3 = 1.5 * log(t);
      if (sv < SV_THRESH) { // no sv
        mult = log_a - LOG_2PI_2 - log_t3 - v * a * w - v * v * t / 2;
        return mult + log(small_sum_2014(t, a, w, ks));
      } else { // sv
        double log_svt = 0.5 * log(1 + sv * t);
        mult = log_a - log_t3 - log_svt + (sv * a * a * w * w
              - 2 * v * a * w - v * v * t) / (2 + 2 * sv * t);
        return mult + log(small_sum_2014(t, a, w, ks));
      }
    } else { // no log
      if (sv < SV_THRESH) { // no sv
        mult = a * exp(-v * a * w - v * v * t / 2) / sqrt(2 * M_PI * t * t * t);
        return mult * small_sum_2014(t, a, w, ks);
      } else { // sv
        mult = a * exp((-v * v * t - 2 * v * a * w + sv * a * a * w * w)
              / (2 + 2 * sv * t)) / sqrt(t * t * t + sv * t * t * t * t);
        return mult * small_sum_2014(t, a, w, ks);
      }
    }
  } else { // large time appx is more efficient
    if (log_prob) { // log
      double log_a2 = 2 * log(a);
      mult =  LOG_PI - log_a2 - v * a * w - v * v * t / 2;
      return mult + log(large_sum_Nav(t, a, w, kl));
    } else { // no log
      mult = M_PI * exp(-v * a * w - v * v * t / 2) / (a * a);
      return mult * large_sum_Nav(t, a, w, kl);
    }
  }
}