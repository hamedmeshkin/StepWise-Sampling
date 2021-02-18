 proc US116a {} {
    global kw pFile deg2Rad deg2Rad2 refList wr kw_twic Energy sec Run kwb kwb_twic Eng_b HelixList phiD   subsec  tiD i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 Prc_362    kwB    Prc_HBond358   Prc_Bond BondToT   wrH rang_bond Kw_Bond  Prc_Ang Prc_Bnd phiBnd phiAng          aid_Ang Alpha     Prc_Bond360 Rc1 Rc2 Rc3 KwD BonD AngD Tor addFoRCe
    global atmcrd ts gFile rang IndexList cmForce  alpha2 rFile tFile table Prc_Helix Prc_Flip  wrB ranG aiD   Prc_113  PI Prc_Linear tor Linear tor step condF condA condT subsec0 tiA  alpha_UpBnd beta_UpBnd alpha_DnBnd beta_DnBnd KwBnd    phiD1 bond358  Kw_Bond360 Ang1 A366 HBond Prc_HBond alpha lDiff aid_Bnd Bond LinAngle addForceAngle   Rc4 Bd Ag  n1 n2 n3 n4
    global Prc_Flip force addForce alpha_Up    beta_Up    alpha_Dn    beta_Dn      oldTS IndF refID angleFixer DelPhi Eng qFile lable  Prc_PhiPsi Prc_Rc hFile addForceA RC phi362   alpha_Dn362    addFORCE stt aid_Ang221B aid_Bnd358   WrA   Ang366 AddForceB angleFixer2 Tor bond BoNd aid_Bnd   aid_Ang BnD Tor366 Tor6 Rc5 Rc6 Tag tiA358 tiA366 tiA359 omega1 omega2 omega3 omega4 omega5 omega6  omega7 omega8 omega9 omega10 omega11 omega12 omegaweight1 omegaweight2 omegaweight3 omegaweight4 omegaweight5 omegaweight6 omegaweight7 omegaweight8 omegaweight9 omegaweight10 omegaweight11 omegaweight12 beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10 beta11 beta12 OmegaWeight359 OmegaWeight366 OmegaWeight358 Alpha116a HBond116a Alpha116b HBond116b Alpha117a HBond117a Alpha117b HBond117b Alpha359 HBond359 Bond359 Ang359  Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355    aid_Ang366 aid_Ang359 aid_Ang358 US116a US359 US366 US358 US116b US117a US117b US366 US355 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b ColVar RCsub AlphaRest EnergyBoundary
    
    set RC  [getdihedral $atmcrd($i1) $atmcrd($i2) $atmcrd($i3) $atmcrd($i4)]
    set RCsub $RC
    
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
    foreach aid_phipsi $Alpha116a {
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
    
    foreach bOnD $HBond116a {
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
    foreach CV $CV116a {
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
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# %% Applying the Biased Potential on the Center of Mass %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
        set DelPhi [angleFixer2 [expr ($RC-$wr)]]
        if {$lable=="FF"} {
            set force  [expr -$kw*$DelPhi*$deg2Rad]
            addForce   $force $i1 $i2 $i3 $i4
        } else {
            if {abs($DelPhi) <= $rang} {
                set force  0.0
                
            } elseif  {$DelPhi > 0.0} {
                set force [expr -$kw*($DelPhi-$rang)*$deg2Rad]
                addForce   $force $i1 $i2 $i3 $i4
                
            } else {
                set force [expr -$kw*($DelPhi+$rang)*$deg2Rad]
                addForce   $force $i1 $i2 $i3 $i4
                
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

proc Exchange116a {oldTag} {
    global kw deg2Rad wr kw_twic Energy sec kwb kwb_twic Eng_b HelixList subsec tiD i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 atmcrd ts rang IndexList aiD tiA HBond alpha Bond LinAngle n1 n2 n3 n4               
    global IndF refID angleFixer Eng lable RC angleFixer2 Tag tiA358 tiA366 tiA359 Alpha116a HBond116a Alpha116b HBond116b Alpha117a HBond117a Alpha117b HBond117b Alpha359 HBond359 Bond359 Ang359  Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355 aid_Ang366 aid_Ang359 aid_Ang358 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b ColVar RCsub Exchange116a

    set DelPhi [angleFixer2 [expr ($RC-$wr)]]
    if {$lable=="FF"} {
        set force  [expr -$kw*$DelPhi*$deg2Rad]
    } else {
        set RC  [getdihedral $atmcrd($i1) $atmcrd($i2) $atmcrd($i3) $atmcrd($i4)]
        set DelPhi [angleFixer2 [expr ($RC-$wr)]]
        if {abs($DelPhi) <= $rang} {
            set force  0.0
        } elseif  {$DelPhi > 0.0} {
            set force [expr -$kw*($DelPhi-$rang)*$deg2Rad]
        } else {
            set force [expr -$kw*($DelPhi+$rang)*$deg2Rad]
        }
    }
    set EnergyA [expr $force * $force / $kw_twic] 
    
   if {$Tag != $oldTag} {
    set Eng_b {}
    foreach aid_phipsi $Alpha116a {
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
    
    
    foreach bOnD $HBond116a {
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
    foreach CV $CV116a {
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

 
