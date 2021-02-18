   proc US117a {} {
    global kw pFile deg2Rad deg2Rad2 refList wr kw_twic Energy sec Run kwb kwb_twic Eng_b HelixList phiD   subsec  tiD i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 Prc_362    kwB    Prc_HBond358   Prc_Bond BondToT   wrH rang_bond Kw_Bond  Prc_Ang Prc_Bnd phiBnd phiAng          aid_Ang Alpha     Prc_Bond360 Rc1 Rc2 Rc3 KwD BonD AngD Tor addFoRCe
    global atmcrd ts gFile rang IndexList cmForce  alpha2 rFile tFile table Prc_Helix Prc_Flip  wrB ranG aiD   Prc_113  PI Prc_Linear tor Linear tor step condF condA condT subsec0 tiA  alpha_UpBnd beta_UpBnd alpha_DnBnd beta_DnBnd KwBnd    phiD1 bond358  Kw_Bond360 Ang1 A366 HBond Prc_HBond alpha lDiff aid_Bnd Bond LinAngle addForceAngle   Rc4 Bd Ag  n1 n2 n3 n4
    global Prc_Flip force addForce alpha_Up    beta_Up    alpha_Dn    beta_Dn      oldTS IndF refID angleFixer DelPhi Eng qFile lable  Prc_PhiPsi Prc_Rc hFile addForceA RC phi362   alpha_Dn362    addFORCE stt aid_Ang221B aid_Bnd358   WrA   Ang366 AddForceB angleFixer2 Tor bond BoNd aid_Bnd   aid_Ang BnD Tor366 Tor6 Rc5 Rc6 Tag tiA358 tiA366 tiA359 omega1 omega2 omega3 omega4 omega5 omega6  omega7 omega8 omega9 omega10 omega11 omega12 omegaweight1 omegaweight2 omegaweight3 omegaweight4 omegaweight5 omegaweight6 omegaweight7 omegaweight8 omegaweight9 omegaweight10 omegaweight11 omegaweight12 beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10 beta11 beta12 OmegaWeight359 OmegaWeight366 OmegaWeight358 Alpha116a HBond116a Alpha116b HBond116b Alpha117a HBond117a Alpha117b HBond117b Alpha359 HBond359 Bond359 Ang359  Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355    aid_Ang366 aid_Ang359 aid_Ang358 US116a US359 US366 US358 US116b US117a US117b US366 US355 tiT358  tiT366 tiT359 tiT8 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b ColVar RCsub aid_Ang117 OmegaWeight117 omega13 omega14 beta13 beta14 beta15 omega15 tiA117 Bound117a AlphaRest EnergyBoundary
    
    set phi1  [expr [getdihedral $atmcrd($i1) $atmcrd($i2) $atmcrd($i3) $atmcrd($i4)] * $deg2Rad]
    set phi2  [expr [getdihedral $atmcrd($j1) $atmcrd($j2) $atmcrd($j3) $atmcrd($j4)] * $deg2Rad]
        
    if {$phi1 > 0.34906585} {set phi1 [expr $phi1 - 6.28318530717959]}
  
    set Rc1 [expr $omega13 * $phi1 - $beta13]
    set Rc2 [expr $omega14 * $phi2 - $beta14]
     

    set RC  [expr 0.5000 * ($Rc1 +  $Rc2)] 
    set RCsub "$phi1 $phi2"
    
    if {$oldTS < $ts} {
        set oldTS $ts
        puts $pFile [format %0.10f $RC]
        if {$ts % 100000 == 0} {
            table
#             puts $rFile "PHiPS =  $Prc_PhiPsi \t\t $Prc_HBond "
            flush $pFile
       }
    }  
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# %% Apply Boundary on Helix resid 361 Phi Psi Angle %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foreach rc_bound $Bound117a subrc $RCsub {T1 T2 T3 T4} $tiA117 {
            lassign $rc_bound  wr_2J rang_RC WrB TwoWrB
            if {$T2!=$T4} {
                set Del_RC   [angleFixer [expr ($subrc-$wr_2J)]]
            } else {
                set Del_RC   [expr ($subrc-$wr_2J)]
            }
            if {abs($Del_RC) <= $rang_RC} {
                set Egg 0.0
            } elseif  {$Del_RC > 0.0} {
                set force_b [expr -$WrB*($Del_RC-$rang_RC)]
                addFoRCe $force_b $subrc $T1 $T2 $T3 $T4
                set Egg [expr $force_b * $force_b / $TwoWrB] 
            } else {
                set force_b [expr -$WrB*($Del_RC+$rang_RC)]
                addFoRCe $force_b $subrc $T1 $T2 $T3 $T4
                set Egg [expr $force_b * $force_b / $TwoWrB]  
            }
    }
    
#     foreach aid_phipsi $AlphaRest {
#             lassign $aid_phipsi T1 T2 T3 T4 wr_2J rang_phipsi WrB 
#             set phiD_phipsi  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
#             set Del_phipsi   [angleFixer2 [expr ($phiD_phipsi-$wr_2J)]]
#             if {abs($Del_phipsi) <= $rang_phipsi} {
#             } elseif  {$Del_phipsi > 0.0} {
#                 set force_b [expr -$WrB*($Del_phipsi-$rang_phipsi)*$deg2Rad]
#                 addFORCE $force_b  $T1 $T2 $T3 $T4
#             } else {
#                 set force_b [expr -$WrB*($Del_phipsi+$rang_phipsi)*$deg2Rad]
#                 addFORCE $force_b  $T1 $T2 $T3 $T4
#             }
#     }



   set Eng_b {}
    foreach aid_phipsi $Alpha117a {
            lassign $aid_phipsi T1 T2 T3 T4 wr_2J rang_phipsi WrB 
            if {$T3!=$T4} {
                set phiD_phipsi  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
            } else {
                set phiD_phipsi  [getangle $atmcrd($T1) $atmcrd($T2) $atmcrd($T3)]
            }
            set Del_phipsi   [angleFixer2 [expr ($phiD_phipsi-$wr_2J)]]
            if {abs($Del_phipsi) <= $rang_phipsi} {
                set Eg 0.0
            } elseif  {$Del_phipsi > 0.0} {
                set force_b [expr -$WrB*($Del_phipsi-$rang_phipsi)*$deg2Rad]
                addFORCE $force_b  $T1 $T2 $T3 $T4
                set Eg   [expr $force_b * $force_b / $WrB / 2.0]
            } else {
                set force_b [expr -$WrB*($Del_phipsi+$rang_phipsi)*$deg2Rad]
                addFORCE $force_b  $T1 $T2 $T3 $T4
                set Eg   [expr $force_b * $force_b / $WrB / 2.0]
            }
            lappend Eng_b $Eg
    }
    
    foreach bOnD $HBond117a {
        lassign $bOnD D1D D2D wrH rang_bond Kw_Bond
        set bond [getbond $atmcrd($D1D) $atmcrd($D2D)]
        set Delbond [expr $bond - $wrH]
        if {abs($Delbond) <= $rang_bond} {
            set Eg 0.0
        } elseif {$Delbond > 0.0} {
            set force_b [expr -$Kw_Bond*($Delbond-$rang_bond)]
            AddForceB [expr $force_b / $bond]  $D1D $D2D
            set Eg [expr $force_b * $force_b / $Kw_Bond / 2.0]
        } else {
            set force_b [expr -$Kw_Bond*($Delbond+$rang_bond)]
            AddForceB [expr $force_b / $bond]  $D1D $D2D
            set Eg [expr $force_b * $force_b / $Kw_Bond / 2.0]
        }
        lappend Eng_b $Eg
    }
    
    set ColVar {}
    foreach CV $CV117a {
            lassign $CV T1 T2 T3 T4 wr_2J rang_cv WrB TwoWrB
            if {$T2!=$T4} {
                set cv  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
                set Del_cv   [angleFixer2 [expr ($cv-$wr_2J)]]
            } else {
                set cv  [getbond $atmcrd($T1) $atmcrd($T2)]
                set Del_cv   [expr ($cv-$wr_2J)]
            }
            lappend ColVar $cv
            if {abs($Del_cv) <= $rang_cv} {
                set Eg 0.0
            } elseif  {$Del_cv > 0.0} {
                set force_b [expr -$WrB*($Del_cv-$rang_cv)]
                addFoRCe $force_b $cv $T1 $T2 $T3 $T4
                set Eg [expr $force_b * $force_b / $TwoWrB] 
            } else {
                set force_b [expr -$WrB*($Del_cv+$rang_cv)]
                addFoRCe $force_b $cv $T1 $T2 $T3 $T4
                set Eg [expr $force_b * $force_b / $TwoWrB]  
            }
            lappend Eng_b $Eg
    }
 
         foreach bnd $aid_Ang117   {
             lassign $bnd  alpha beta gamma Delta KwL
             set AB [expr  $alpha * exp($beta  * $RC)]
             set CD [expr  $gamma * exp($Delta * $RC)]
             set Y  [expr $AB+$CD]
             set const [expr  $beta * $AB + $Delta * $CD]
             set coefP [expr -2.0 + $const]
             set coefM [expr  2.0 + $const]
         }
         
         set DelBond  [expr $Rc1 - $Rc2]
         if {abs($DelBond) <= $Y} {
             set Eg 0.0
         } elseif  {$DelBond > 0.0} {
             set ForceAngle   [expr  $KwL * ($DelBond - $Y)]
             set Fe  [expr $ForceAngle * $const]
             set FfP [expr $ForceAngle * $coefP]
             set FfM [expr $ForceAngle * $coefM]
             set Forces "$FfP $FfM"
             foreach omegaweight $OmegaWeight117 subrc $RCsub FM $Forces {T1 T2 T3 T4} $tiA117 {
                 addFoRCe [expr $omegaweight * $FM] $subrc $T1 $T2 $T3 $T4
             }
         } else {
             set ForceAngle   [expr -$KwL * ($DelBond + $Y)]
             set Fe  [expr $ForceAngle * $const]
             set FfM [expr $ForceAngle * $coefM]
             set FfP [expr $ForceAngle * $coefP]
             set Forces "$FfM $FfP"
             foreach omegaweight $OmegaWeight117 subrc $RCsub FM $Forces {T1 T2 T3 T4} $tiA117 {
                 addFoRCe [expr $omegaweight * $FM] $subrc $T1 $T2 $T3 $T4
             }
         }
        
 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# %% Applying the Biased Potential on the Center of Mass %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
        set DelPhi [expr ($RC-$wr)] 
        if {$lable=="FF"} {
            set force  [expr  -$kw*$DelPhi]
            foreach omegaweight $OmegaWeight117 subrc $RCsub {T1 T2 T3 T4} $tiA117 {
                addFoRCe   [expr $omegaweight * $force] $subrc $T1 $T2 $T3 $T4
            }
        } else {
            if {abs($DelPhi) <= $rang} {
                set force  0.0
                set Eng    0.0
            } elseif  {$DelPhi > 0.0} {
                set force [expr -$kw*($DelPhi-$rang)]
                foreach omegaweight $OmegaWeight117 subrc $RCsub {T1 T2 T3 T4} $tiA117 {
                    addFoRCe   [expr $omegaweight * $force] $subrc $T1 $T2 $T3 $T4
                }
            } else {
                set force [expr -$kw*($DelPhi+$rang)]
                foreach omegaweight $OmegaWeight117 subrc $RCsub {T1 T2 T3 T4} $tiA117 {
                    addFoRCe   [expr $omegaweight * $force] $subrc $T1 $T2 $T3 $T4
                }
            }
        }
        set Energy [expr $force * $force / $kw_twic] 
    
    if {$ts % 200 == 0} {    
        if {$lable == "FA"} {
            set EnergyBoundary [expr [join $Eng_b +]]
            set Energy [expr $Energy + $EnergyBoundary]
        }
    }
return
}


proc Exchange117a {oldTag} {
    global kw deg2Rad wr kw_twic Energy sec kwb kwb_twic Eng_b HelixList subsec tiD i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 atmcrd ts rang IndexList aiD tiA HBond alpha Bond LinAngle n1 n2 n3 n4               
    global IndF refID angleFixer Eng lable RC angleFixer2 Tag tiA358 tiA366 tiA359 Alpha116a HBond116a Alpha116b HBond116b Alpha117a HBond117a Alpha117b HBond117b Alpha359 HBond359 Bond359 Ang359  Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355 aid_Ang366 aid_Ang359 aid_Ang358 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b ColVar RCsub Exchange116a Exchange117a OmegaWeight117 omega13 omega14 beta13 beta14 beta15 omega15
    
    
    
        if {$lable=="FF"} {
            set DelPhi [expr ($RC-$wr)] 
            set force  [expr  -$kw*$DelPhi]
        } else {
            set phi1  [expr [getdihedral $atmcrd($i1) $atmcrd($i2) $atmcrd($i3) $atmcrd($i4)] * $deg2Rad]
            set phi2  [expr [getdihedral $atmcrd($j1) $atmcrd($j2) $atmcrd($j3) $atmcrd($j4)] * $deg2Rad]
            if {$phi1 > 0.34906585} {set phi1 [expr $phi1 - 6.28318530717959]}
            set Rc1 [expr $omega13 * $phi1 - $beta13]
            set Rc2 [expr $omega14 * $phi2 - $beta14]
            set RC  [expr 0.500 * ($Rc1 +  $Rc2)] 
            set DelPhi [expr ($RC-$wr)] 
            if {abs($DelPhi) <= $rang} {
                set force  0.0
                set Eng    0.0
            } elseif  {$DelPhi > 0.0} {
                set force [expr -$kw*($DelPhi-$rang)]
            } else {
                set force [expr -$kw*($DelPhi+$rang)]
            }
        }
       set EnergyA [expr $force * $force / $kw_twic] 

    if {$Tag != $oldTag} {
    set Eng_b {}
    foreach aid_phipsi $Alpha117a {
            lassign $aid_phipsi T1 T2 T3 T4 wr_2J rang_phipsi WrB 
            if {$T3!=$T4} {
                set phiD_phipsi  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
            } else {
                set phiD_phipsi  [getangle $atmcrd($T1) $atmcrd($T2) $atmcrd($T3)]
            }
            set Del_phipsi   [angleFixer2 [expr ($phiD_phipsi-$wr_2J)]]
            if {abs($Del_phipsi) <= $rang_phipsi} {
                set Eg 0.0
            } elseif  {$Del_phipsi > 0.0} {
                set DEl [expr ($Del_phipsi-$rang_phipsi)*$deg2Rad]
                set Eg  [expr $DEl * $DEl * $WrB / 2.0]
            } else {
                set DEl [expr ($Del_phipsi-$rang_phipsi)*$deg2Rad]
                set Eg  [expr $DEl * $DEl * $WrB / 2.0]
            }
            lappend Eng_b $Eg
    }
    
    foreach bOnD $HBond117a {
        lassign $bOnD D1D D2D wrH rang_bond Kw_Bond
        set bond [getbond $atmcrd($D1D) $atmcrd($D2D)]
        set Delbond [expr $bond - $wrH]
        if {abs($Delbond) <= $rang_bond} {
            set Eg 0.0
        } elseif {$Delbond > 0.0} {
            set DEl [expr ($Delbond-$rang_bond)]
            set Eg  [expr $DEl * $DEl * $Kw_Bond / 2.0]
        } else {
            set DEl [expr ($Delbond-$rang_bond)]
            set Eg  [expr $DEl * $DEl * $Kw_Bond / 2.0]
        }
        lappend Eng_b $Eg
    }
    
    set ColVar {}
    foreach CV $CV117a {
            lassign $CV T1 T2 T3 T4 wr_2J rang_cv WrB TwoWrB
            if {$T2!=$T4} {
                set cv  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
                set Del_cv   [angleFixer2 [expr ($cv-$wr_2J)]]
            } else {
                set cv  [getbond $atmcrd($T1) $atmcrd($T2)]
                set Del_cv   [expr ($cv-$wr_2J)]
            }
            lappend ColVar $cv
            if {abs($Del_cv) <= $rang_cv} {
                set Eg 0.0
            } elseif  {$Del_cv > 0.0} {
                set force_b [expr -$WrB*($Del_cv-$rang_cv)]
                set Eg [expr $force_b * $force_b / $TwoWrB] 
            } else {
                set force_b [expr -$WrB*($Del_cv+$rang_cv)]
                set Eg [expr $force_b * $force_b / $TwoWrB]  
            }
            lappend Eng_b $Eg
    }
   
    set EnergyA [expr $EnergyA + [expr [join $Eng_b +]]]
   }
    
    return $EnergyA
}

 
