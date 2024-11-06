// -*-Faust-*-

declare name "leveler";
declare version "0.1";
declare author "Bart Brouns";
declare license "GPLv3";

// double precision -double needed!

ebu = library("ebur128.lib");
ex = library("expanders.lib");
import("stdfaust.lib");

process = leveler_sc(target)~(_, _);

///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                    //
///////////////////////////////////////////////////////////////////////////////

leveler_meter_gain = hbargraph("v:/[1][unit:dB]gain",-50,50);
bp = checkbox("v:/[2]bypass") : si.smoo;
target = hslider("v:/[3]target[unit:dB]", init_leveler_target,-50,-2,1);
leveler_speed = hslider("v:/[4][unit:%][integer]speed", init_leveler_speed, 0, 100, 1) * 0.0015; //.005, 0.15, .005);
leveler_brake_thresh = /*target + */hslider("v:/[5][unit:dB]brake threshold", init_leveler_brake_threshold,-90,0,1);
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

lk2_var(Tg)= par(i,2,kfilter : zi) :> 4.342944819 * log(max(1e-12)) : -(0.691) with {
  // maximum assumed sample rate is 192k
  maxSR = 192000;
  sump(n) = ba.slidingSump(n, Tg*maxSR)/max(n,ma.EPSILON);
  envelope(period, x) = x * x :  sump(rint(period * ma.SR));
  zi = envelope(Tg); // mean square: average power = energy/Tg = integral of squared signal / Tg

  kfilter = ebu.prefilter;
};

lk2 = lk2_var(3);
lk2_short = lk2_var(0.4);
lufs_meter(l,r) = l,r <: l, attach(r, (lk2 : vbargraph("[unit:dB]out-lufs-s",-120,0))) : _,_;

///////////////////////////////////////////////////////////////////////////////
//                                  LEVELER                                  //
///////////////////////////////////////////////////////////////////////////////

leveler_sc(target,fl,fr,l,r) =
  (calc(lk2_short(fl,fr))*(1-bp)+bp)
  <: (_*l,_*r)
with {
  lp1p(cf) = si.smooth(ba.tau2pole(1/(2*ma.PI*cf)));
  calc(lufs) = FB(lufs)~_: ba.db2linear;
  FB(lufs,prev_gain) =
    (target - lufs)
    +(prev_gain )
    :  limit(limit_neg,limit_pos)
    : lp1p(leveler_speed_brake(abs(l)+abs(r)))
    : leveler_meter_gain;

  limit(lo,hi) = min(hi) : max(lo);

  leveler_speed_brake(sc) = (expander(sc) <: attach(_, (1-_) : meter_leveler_brake)) : _ * leveler_speed;

  expander(x) = (ex.peak_expansion_gain_mono_db(maxHold,strength,leveler_brake_thresh,range,gate_att,hold,gate_rel,knee,prePost,x)
                 : ba.db2linear
                 :max(0)
                 :min(1));

  maxHold = hold*192000;
  strength = 2;
  range = -120;
  gate_att = 0.05;
  hold = 0.1;
  gate_rel = 0.3;
  knee = 12;
  prePost = 1;
};
