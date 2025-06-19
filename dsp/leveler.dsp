// -*-Faust-*-

declare name "leveler";
declare version "0.1";
declare author "Bart Brouns";
declare license "GPLv3";

// double precision -double needed!

ebu = library("ebur128.lib");
ex = library("expanders.lib");
ds = library("dynamic-smoothing.lib");
import("stdfaust.lib");

maxSR = 192000;
// maxSR = 48000;

process(l,r) =

  ( ((l,r):leveler_sc(target)~(_,_)
                              :(
       (_*(1-bp))
      ,(_*(1-bp))
     ))
  , (l*bp,r*bp)
  ):>(_,_);

///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                    //
///////////////////////////////////////////////////////////////////////////////

basefreq =
  it.interpolate_linear(leveler_speed
                        :pow(
                          2 // hslider("base freq power", 2, 0.1, 10, 0.1)
                        )
                       , 0.01
                       , 0.2 // hslider("base freq fast", 0.2, 0.1, 0.3, 0.001)
                       );

sensitivity =
  it.interpolate_linear(leveler_speed
                        :pow(
                          0.5 // hslider("sens power", 0.5, 0.1, 10, 0.1)
                        )
                       , 0.00000025
                       , 0.0000025 // hslider("sens fast", 0.0000025, 0.0000025, 0.000005, 0.0000001)
                       );
cf = hslider("smoo freq", 0.1, 0.01, 1000, 0.01);
low_freq = hslider("low_freq", 0.13, 0.01, 10, 0.01);
hi_freq = hslider("hi_freq", 2.5, 1, 10000, 0.5);

leveler_meter_gain = hbargraph("v:/[1][unit:dB]gain",-50,50);
bp = checkbox("v:/[2]bypass"):si.smoo;
target = hslider("v:/[3]target[unit:dB]", init_leveler_target,-50,-2,1);
leveler_speed = hslider("v:/[4][unit:%][integer]speed", init_leveler_speed, 0, 100, 1) * 0.01; //.005, 0.15, .005);
leveler_brake_thresh = target + hslider("v:/[5][unit:dB]brake threshold", init_leveler_brake_threshold,-90,0,1)+32;
meter_leveler_brake = _*100 : hbargraph("v:/[6][unit:%][integer]brake",0,100);
limit_pos = hslider("v:/[7][unit:dB]max boost", init_leveler_maxboost, 0, 60, 1);
limit_neg = hslider("v:/[8][unit:dB]max cut", init_leveler_maxcut, 0, 60, 1) : ma.neg;

init_leveler_target = -18;
init_leveler_maxboost = 20;
init_leveler_maxcut = 20;
init_leveler_brake_threshold = -14;
init_leveler_speed = 20;

///////////////////////////////////////////////////////////////////////////////
//                                 LUFS METER                                //
///////////////////////////////////////////////////////////////////////////////

lk2_fixed(Tg)= par(i,2,kfilter : zi) :> 4.342944819 * log(max(1e-12)) : -(0.691) with {
  // maximum assumed sample rate is 192k
  sump(n) = ba.slidingSump(n, Tg*maxSR)/max(n,ma.EPSILON);
  envelope(period, x) = x * x :  sump(rint(period * ma.SR));
  zi = envelope(Tg); // mean square: average power = energy/Tg = integral of squared signal / Tg

  kfilter = ebu.prefilter;
};

lk2_var(Tg)= par(i,2,kfilter : zi) :> 4.342944819 * log(max(1e-12)) : -(0.691) with {
  // maximum assumed sample rate is 192k
  sump(n) = ba.slidingSump(n, 0.4 * maxSR)/max(n,ma.EPSILON);
  envelope(period, x) = x * x :  sump(rint(period * ma.SR));
  zi = envelope(Tg); // mean square: average power = energy/Tg = integral of squared signal / Tg

  kfilter = ebu.prefilter;
};
lk2 = lk2_fixed(3);
lk2_short = lk2_fixed(0.4);
lufs_meter(l,r) = l,r <: l, attach(r, (lk2 : vbargraph("[unit:dB]out-lufs-s",-120,0))) : _,_;

///////////////////////////////////////////////////////////////////////////////
//                                  LEVELER                                  //
///////////////////////////////////////////////////////////////////////////////

lk2_time =
  // 0.4;
  hslider("lk2 time", 0.01, 0.001, 3, 0.001);
// it.interpolate_linear(leveler_speed :pow(hslider("lk2 power", 2, 0.1, 10, 0.1))
//                      ,0.4 // hslider("lk2 time", 0.4, 0.001, 3, 0.001)
//                      , 0.04):max(0);
leveler_sc(target,fl,fr,l,r) =
  calc(lk2_fixed(0.01,fl,fr))
  // (calc(lk2_var(lk2_time,fl,fr))*(1-bp)+bp)
  <: (_*l,_*r)
with {
  // lp1p(cf) = si.smooth(ba.tau2pole(1/(2*ma.PI*cf)));
  calc(lufs) = FB(lufs)~_: ba.db2linear;
  FB(lufs,prev_gain) =
    (target - lufs)
    +(prev_gain )
    : ds.dynamicSmoothing(
      sensitivity * expander(abs(fl)+abs(fr))
    ,  basefreq * expander(abs(fl)+abs(fr))
    )
    :  limit(limit_neg,limit_pos)
    : leveler_meter_gain;

  limit(lo,hi) = min(hi) : max(lo);

  leveler_speed_brake(sc) = expander(sc) * leveler_speed;

  // Brake idea and implementation: Klaus Scheuermann⁩,  https://github.com/trummerschlunk

  expander(x) = (ex.peak_expansion_gain_mono_db(maxHold,strength,leveler_brake_thresh,range,gate_att,hold,gate_rel,knee,prePost,x)
                 : ba.db2linear
                 :max(0)
                 :min(1))
                <: attach(_, (1-_) : meter_leveler_brake) ;

  maxHold = hold*maxSR;
  strength = 1;
  // hslider("gate strength", 1, 0.1, 10, 0.1);
  range = 0-(ma.MAX);
  gate_att =
    0;
  // hslider("gate att", 0.0, 0.0, 1, 0.001);
  hold = 0.0001;
  gate_rel =
    0.1;
  // hslider("gate rel", 0.1, 0.0, 1, 0.001);
  knee =
    ma.EPSILON;
  // hslider("gate knee", 0, 0, 90, 1);
  prePost = 1;
};
