
# ################################################################ 
# ##### Reading From Files as initial parameter and so on ######## 
# ################################################################
proc table {} {
    global IndexList refList HelixList Linear refGroup HBond Ang366 Alpha aid_Ang  Tor366 omega1 omega2 omega3 omega4 omega5 omega6  omega7 omega8 omega9 omega10 omega11 omega12 omegaweight1 omegaweight2 omegaweight3 omegaweight4 omegaweight5 omegaweight6 omegaweight7 omegaweight8 omegaweight9 omegaweight10 omegaweight11 omegaweight12  beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10 beta11 beta12 OmegaWeight359 OmegaWeight366 OmegaWeight358 Alpha116a HBond116a Alpha117a HBond117a Alpha117b HBond117b Alpha116b HBond116b Alpha359 HBond359 Bond359 Ang359 Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355 aid_Ang366 aid_Ang359 aid_Ang358 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b aid_Ang117 OmegaWeight117 omega13 omega14 beta13 beta14 beta15 omega15 Bound117a Bound117b Bound359 Bound358 Bound366 Bond358 AlphaRest OmegaWeight117b omega15 omega16 omega17 beta15 beta16 beta17
    
    
    
    set IndexList {}
    set fd [open ../input/torsion_index.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 4} {
            error "Every line must have 4 elements!"
        }
        lappend IndexList $line
        foreach aid $line {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Bound359 {}
    set fd [open ../input/Bound359.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
        lappend Bound359 $line
        lassign $line ind1 ind2 ind3 ind4
    }
    close $fd
    
    set Bound117a {}
    set fd [open ../input/Bound117a.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
        lappend Bound117a $line
        lassign $line ind1 ind2 ind3 ind4
    }
    close $fd
    
    set Bound117b {}
    set fd [open ../input/Bound117b.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
        lappend Bound117b $line
        lassign $line ind1 ind2 ind3 ind4
    }
    close $fd
    
    set Bound366 {}
    set fd [open ../input/Bound366.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
        lappend Bound366 $line
        lassign $line ind1 ind2 ind3 ind4
    }
    close $fd
    
    set Bound358 {}
    set fd [open ../input/Bound358.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
        lappend Bound358 $line
        lassign $line ind1 ind2 ind3 ind4
    }
    close $fd
    
    set CV366 {}
    set fd [open ../input/CV366.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV366 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV358 {}
    set fd [open ../input/CV358.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV358 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV359 {}
    set fd [open ../input/CV359.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV359 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV355 {}
    set fd [open ../input/CV355.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV355 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV117b {}
    set fd [open ../input/CV117b.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV117b $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV117a {}
    set fd [open ../input/CV117a.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV117a $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV116b {}
    set fd [open ../input/CV116b.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV116b $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set CV116a {}
    set fd [open ../input/CV116a.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 9} {
            error "Every line must have 9 elements!"
        }
        lappend CV116a $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
  
    set AlphaRest {}
    set fd [open ../input/AlphaRest.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend AlphaRest $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha116a {}
    set fd [open ../input/PHI_PSI354to357_Alpha116a.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha116a $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha116b {}
    set fd [open ../input/PHI_PSI354to357_Alpha116b.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha116b $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha117a {}
    set fd [open ../input/PHI_PSI354to357_Alpha117a.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha117a $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha117b {}
    set fd [open ../input/PHI_PSI354to357_Alpha117b.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha117b $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha359 {}
    set fd [open ../input/PHI_PSI354to357_Alpha359.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha359 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha366 {}
    set fd [open ../input/PHI_PSI354to357_Alpha366.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha366 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha358 {}
    set fd [open ../input/PHI_PSI354to357_Alpha358.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha358 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd
    
    set Alpha355 {}
    set fd [open ../input/PHI_PSI354to357_Alpha355.dat r]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 8} {
            error "Every line must have 8 elements!"
        }
        lappend Alpha355 $line
        lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
        set lin "$ind1 $ind2 $ind3 $ind4"
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    }
    close $fd

    set HBond116a {} 
    set fd [open ../input/HBond116a.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond116a "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond116b {} 
    set fd [open ../input/HBond116b.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond116b "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond117a {} 
    set fd [open ../input/HBond117a.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond117a "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond117b {} 
    set fd [open ../input/HBond117b.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond117b "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond359 {} 
    set fd [open ../input/HBond359.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond359 "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond366 {} 
    set fd [open ../input/HBond366.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond366 "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond358 {} 
    set fd [open ../input/HBond358.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond358 "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd
    
    set HBond355 {} 
    set fd [open ../input/HBond355.dat]
    set OK 1
    while {1} {
        if {[gets $fd line] == -1} {
            break
        }
        if {[llength $line] != 5} {
            error "Every line must have 5 elements!"
        }
    
        lassign $line ind1 ind2 ind3 ind4 ind5 
        lappend HBond355 "$ind1 $ind2 $ind3 $ind4 $ind5"
        set lin "$ind1 $ind2" 
        foreach aid $lin {
            if {![catch {expr {abs($aid)}}]} {
                if {![info exists flag($aid)]} {
                    addatom $aid
                    set flag($aid) 1
                }
            }
        }
    } 
    close $fd

set Ang359 {} 
set fd [open ../input/Ang359.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 8} {
        error "Every line must have 8 elements!"
    }
    lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8   
    lappend Ang359 "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6 $ind7 $ind8"
    set lin "$ind1 $ind2" 
    foreach aid $lin {
        if {![catch {expr {abs($aid)}}]} {
            if {![info exists flag($aid)]} {
                addatom $aid
                set flag($aid) 1
            }
        }
    }
}
close $fd

set Bond358 {} 
set fd [open ../input/Bond358.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 6} {
        error "Every line must have 6 elements!"
    }
    lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
    lappend Bond358 "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6"
    set lin "$ind1 $ind2" 
    foreach aid $lin {
        if {![catch {expr {abs($aid)}}]} {
            if {![info exists flag($aid)]} {
                addatom $aid
                set flag($aid) 1
            }
        }
    }
}
close $fd  
    
set Bond359 {} 
set fd [open ../input/Bond359.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 6} {
        error "Every line must have 6 elements!"
    }
    lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
    lappend Bond359 "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6"
    set lin "$ind1 $ind2" 
    foreach aid $lin {
        if {![catch {expr {abs($aid)}}]} {
            if {![info exists flag($aid)]} {
                addatom $aid
                set flag($aid) 1
            }
        }
    }
}
close $fd   

set Bond366 {} 
set fd [open ../input/Bond366.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 6} {
        error "Every line must have 6 elements!"
    }
    lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8
    lappend Bond366 "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6"
    set lin "$ind1 $ind2" 
    foreach aid $lin {
        if {![catch {expr {abs($aid)}}]} {
            if {![info exists flag($aid)]} {
                addatom $aid
                set flag($aid) 1
            }
        }
    }
}
close $fd 
 
set aid_Ang359 {} 
set fd [open ../input/aid_Ang359.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 5} {
        error "Every line must have 5 elements!"
    }

    lassign $line ind1 ind2 ind3 ind4 ind5
    lappend aid_Ang359 "$ind1 $ind2 $ind3 $ind4 $ind5" 
} 
close $fd

set aid_Ang358 {} 
set fd [open ../input/aid_Ang358.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 5} {
        error "Every line must have 5 elements!"
    }

    lassign $line ind1 ind2 ind3 ind4 ind5
    lappend aid_Ang358 "$ind1 $ind2 $ind3 $ind4 $ind5" 
} 
close $fd

set aid_Ang366 {} 
set fd [open ../input/aid_Ang366.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 11} {
        error "Every line must have 11 elements!"
    }

    lassign $line  ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8 ind9 ind10 ind11
    lappend aid_Ang366 "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6 $ind7 $ind8 $ind9 $ind10 $ind11" 
} 
close $fd

set aid_Ang117 {} 
set fd [open ../input/aid_Ang117.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 5} {
        error "Every line must have 5 elements!"
    }

    lassign $line ind1 ind2 ind3 ind4 ind5
    lappend aid_Ang117 "$ind1 $ind2 $ind3 $ind4 $ind5" 
} 
close $fd

set IniFun {}
set fd [open ../input/Func.dat]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 3} {
        error "Every line must have 3 elements!"
    }
    lappend IniFun $line   
}
close $fd

lassign [lindex $IniFun 0] w1 Delphi1 phiI1
lassign [lindex $IniFun 1] w2 Delphi2 phiI2
lassign [lindex $IniFun 2] w3 Delphi3 phiI3
lassign [lindex $IniFun 3] w4 Delphi4 phiI4
lassign [lindex $IniFun 4] w5 Delphi5 phiI5
lassign [lindex $IniFun 5] w6 Delphi6 phiI6
lassign [lindex $IniFun 6] w7 Delphi7 phiI7
lassign [lindex $IniFun 7] w8 Delphi8 phiI8
lassign [lindex $IniFun 8] w9 Delphi9 phiI9
lassign [lindex $IniFun 9] w10 Delphi10 phiI10
lassign [lindex $IniFun 10] w11 Delphi11 phiI11
lassign [lindex $IniFun 11] w12 Delphi12 phiI12
lassign [lindex $IniFun 12] w13 Delphi13 phiI13
lassign [lindex $IniFun 13] w14 Delphi14 phiI14
lassign [lindex $IniFun 14] w15 Delphi15 phiI15
lassign [lindex $IniFun 15] w16 Delphi16 phiI16
lassign [lindex $IniFun 16] w17 Delphi17 phiI17

foreach "omega1 omega2 omega3" "[expr 1.0/$Delphi1] [expr 1.0/$Delphi2] [expr 1.0/$Delphi3]" {}
foreach "omega4 omega5 omega6  omega7 omega8" "[expr 1.0/$Delphi4] [expr 1.0/$Delphi5] [expr 1.0/$Delphi6] [expr 1.0/$Delphi7] [expr 1.0/$Delphi8]" {}
foreach "omega9 omega10 omega11 omega12" "[expr 1.0/$Delphi9] [expr 1.0/$Delphi10] [expr 1.0/$Delphi11] [expr 1.0/$Delphi12]" {}
foreach "omega13 omega14" "[expr 1.0/$Delphi13] [expr 1.0/$Delphi14]" {}
foreach "omega15 omega16 omega17"  "[expr 1.0/$Delphi15] [expr 1.0/$Delphi16] [expr 1.0/$Delphi17]" {}

foreach "omegaweight1 omegaweight2 omegaweight3" "[expr $w1/$Delphi1] [expr $w2/$Delphi2] [expr $w3/$Delphi3]" {}
foreach "omegaweight4 omegaweight5 omegaweight6 omegaweight7 omegaweight8" "[expr $w4/$Delphi4] [expr $w5/$Delphi5] [expr $w6/$Delphi6] [expr $w7/$Delphi7] [expr $w8/$Delphi8]" {}
foreach "omegaweight9 omegaweight10 omegaweight11 omegaweight12" "[expr $w9/$Delphi9] [expr $w10/$Delphi10] [expr $w11/$Delphi11] [expr $w12/$Delphi12]" {}
foreach "omegaweight13 omegaweight14 omegaweight15" "[expr $w13/$Delphi13] [expr $w14/$Delphi14] [expr $w15/$Delphi15]" {}
foreach "omegaweight15 omegaweight16 omegaweight17" "[expr $w15/$Delphi15] [expr $w16/$Delphi16] [expr $w17/$Delphi17]" {}

foreach "beta1 beta2 beta3" "[expr $omega1*$phiI1] [expr $omega2*$phiI2] [expr $omega3*$phiI3]" {}
foreach "beta4 beta5 beta6 beta7 beta8" "[expr $omega4*$phiI4] [expr $omega5*$phiI5] [expr $omega6*$phiI6] [expr $omega7*$phiI7] [expr $omega8*$phiI8]" {}
foreach "beta9 beta10 beta11 beta12" "[expr $omega9*$phiI9] [expr $omega10*$phiI10] [expr $omega11*$phiI11] [expr $omega12*$phiI12]" {}
foreach "beta13 beta14" "[expr $omega13*$phiI13] [expr $omega14*$phiI14]" {}
foreach "beta15 beta16 beta17" "[expr $omega15*$phiI15] [expr $omega16*$phiI16] [expr $omega17*$phiI17]" {}

set OmegaWeight359 " $omegaweight2  $omegaweight3"
set OmegaWeight366 "$omegaweight4  $omegaweight5  $omegaweight6  $omegaweight7 $omegaweight8"
set OmegaWeight358 "$omegaweight9  $omegaweight10 $omegaweight11"
set OmegaWeight117 "$omegaweight13 $omegaweight14"
set OmegaWeight117b "$omegaweight15 $omegaweight16 $omegaweight17"

puts $IniFun
 
    set fd [open ../input/ref_win0.dat r]
    set refList [read $fd]
    close $fd
    
    set fd [open ../input/Group.dat r]
    set refGroup [read $fd]
    close $fd


    return
}
##--------------------------------------------------------------------------##
table

set groupid1  [addgroup [lindex $refGroup 0]]
set groupid2  [addgroup [lindex $refGroup 1]]
set groupid3  [addgroup [lindex $refGroup 2]]
set groupid4  [addgroup [lindex $refGroup 3]]
set groupid5  [addgroup [lindex $refGroup 4]]
set groupid6  [addgroup [lindex $refGroup 5]]
set groupid7  [addgroup [lindex $refGroup 6]]
set groupid8  [addgroup [lindex $refGroup 7]]
set groupid9  [addgroup [lindex $refGroup 8]]
set groupid10 [addgroup [lindex $refGroup 9]]
set groupid11 [addgroup [lindex $refGroup 10]]
set groupid12 [addgroup [lindex $refGroup 11]]

set nRef [llength $refList]
# ############################################################################# 
# ##### Testin the Metropolis Criterion to Calculate Energy Difference ######## 
# #############################################################################
proc calcDiffEng {n} {  
    global refID IndF angleFixer deg2Rad2 atmcrd nRef wr lable HelixList tFile kw sec DelPhi ts qFile rang Energy rFile deg2Rad kw_twic kwb_twic kwb phiD1 wrB ranG aiD Eng subsec0 i1 i2 i3 i4 GetDihedAngle RC  omega1 omega2 omega3 Rc359 Rc362 Rc363 Rc366 Rc367 phi113 BondToT phiD tor condF condA condT step aid_Bnd360 phiBnd phiAng Ang221B  Tor430  Ang221S  bond221   phiD1 bond358 phiAngTot A366 BoNd tor lDiff Bond Rc1 Rc2 Rc3 BonD AngD Tor Ltor BnD Tor6 Rc4 Rc5 Rc6 Bd Ag RCsub ColVar Exchange116a Exchange116b Exchange117a  Exchange117b Exchange355 Exchange358 Exchange359 Exchange366 Tag OmegaWeight117 omega13 omega14 beta13 beta14 beta15 omega15 EnergyBoundary Rcm
    
    puts $qFile "$ts $RC $sec $kw $RCsub $ColVar"
    # If n is even, window 2*i exchanges with 2*i + 1;
    # If n is odd, window 2*i exchanges with 2*i - 1.
    if {($refID + $n) % 2 == 0} {
        set newRef [expr $refID + 1]
        if {$newRef >= $nRef} {
            return nan 
        }
    } else {
        set newRef [expr $refID - 1]
        if {$newRef < 0} {
            return nan
        }
    }
 
    set oldTag $Tag 
    set oldlable $lable
    IndF $newRef
    
    if {$oldlable != $lable} {
        if {$oldlable == "FA"} {
            set Energy [expr $Energy - $EnergyBoundary]
        }
    }
    
    if {$Tag == "R117a"} {
        return [expr [Exchange117a $oldTag]-$Energy]  
    } elseif {$Tag == "R116a"} {
        return [expr [Exchange116a $oldTag]-$Energy] 
    } elseif {$Tag == "R359"}  {
        return [expr [Exchange359 $oldTag]-$Energy] 
    } elseif {$Tag == "R116b"} {
        return [expr [Exchange116b $oldTag]-$Energy] 
    } elseif {$Tag == "R366"}  {
        return [expr [Exchange366 $oldTag]-$Energy] 
    } elseif {$Tag == "R117b"} {
        return [expr [Exchange117b $oldTag]-$Energy] 
    } elseif {$Tag == "R358"}  {
        return [expr [Exchange358 $oldTag]-$Energy] 
    } elseif {$Tag == "R355"}  {
        return [expr [Exchange355 $oldTag]-$Energy]  
    }
 
} 
##-------------------------------------------------------------------##

# ##############################################################################
# ##### Index of the Bias Potential and Assigning the related parameter ######## 
# ##############################################################################
proc IndF {refID} {
    global rang refList lable IndexList kw sec wr condF condA condT kw_twic kwb_twic kwb i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 wrB ranG aiD subsec tiD S1 S2 S3 S4  groupid2 groupid3  R1 R2 R3 R4 step   tiA subsec0 aid_Bnd360 tiG tiQ deg2Rad n1 n2 n3 n4 Tag tiA358 tiA366 tiA359 tiA117 tiA117b
   
    set sec        $refID;
     
    set lable      [lindex $refList $sec 0]
    set Tag        [lindex $refList $sec 1]
    set subsec0    [lindex $refList $sec 2]
    set subsec1    [lindex $refList $sec 3]
    set subsec2    [lindex $refList $sec 4]
    set subsec3    [lindex $refList $sec 5]
    set subsec4    [lindex $refList $sec 6]
   
    set wr         [lindex $refList $sec 7]
    set rang       [lindex $refList $sec 8]
    set kw         [lindex $refList $sec 9]
    
    set aid0     [lindex $IndexList  $subsec0] 
    lassign $aid0 i1 i2 i3 i4
    set aid1     [lindex $IndexList  $subsec1] 
    lassign $aid1 j1 j2 j3 j4
    set aid2     [lindex $IndexList  $subsec2] 
    lassign $aid2 k1 k2 k3 k4
    set aid3     [lindex $IndexList  $subsec3] 
    lassign $aid3 l1 l2 l3 l4
    set aid4     [lindex $IndexList  $subsec4] 
    lassign $aid4 m1 m2 m3 m4
    
    set tiA117 "$aid0 $aid1"
    set tiA117b "$aid0 $aid1 $aid2"
    set tiA359 "$aid1 $aid2"
    set tiA366 "$aid0 $aid1"
    set tiA358 "$aid0 $aid1 $aid2"
    set kw_twic    [expr  $kw*2.0]
}
##------------------------------------------------------------##

# ##################################################################### 
# ##### Fix the Delta phi not to get a value bigger that |180| ######## 
# ##################################################################### 
proc angleFixer {DelPhi} {
    if {$DelPhi > 3.14159265358979} {
        set DelPhi [expr $DelPhi-6.28318530717959]
    } elseif {$DelPhi < -3.14159265358979} {
        set DelPhi [expr $DelPhi+6.28318530717959]
    }
    return $DelPhi
}

proc angleFixer2 {DelPhi} {
    if {$DelPhi > 180} {
        set DelPhi [expr $DelPhi-360]
    } elseif {$DelPhi < -180} {
        set DelPhi [expr $DelPhi+360]
    }
    return $DelPhi
}
##----------------------------------------------------------------------------##

# ############################################### 
# ##### Angle & Aihedral Force On Torsion ####### 
# ###############################################
proc addForce {force T1 T2 T3 T4} {
    global atmcrd ts
    foreach {g1 g2 g3 g4} [dihedralgrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)] {}
    addforce $T1 [vecscale $g1 $force]
    addforce $T2 [vecscale $g2 $force]
    addforce $T3 [vecscale $g3 $force]
    addforce $T4 [vecscale $g4 $force]   
} 
##------------------------------------------------------------##

proc addFORCE {force T1 T2 T3 T4} {
    global atmcrd ts
    if {$T3!=$T4} {
        foreach {g1 g2 g3 g4} [dihedralgrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)] {}
        addforce $T1 [vecscale $g1 $force]
        addforce $T2 [vecscale $g2 $force]
        addforce $T3 [vecscale $g3 $force]
        addforce $T4 [vecscale $g4 $force] 
    } else {
        foreach {g1 g2 g3} [anglegrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3)] {}
        addforce $T1 [vecscale $g1 $force]
        addforce $T2 [vecscale $g2 $force]
        addforce $T3 [vecscale $g3 $force]
    }
} 
##------------------------------------------------------------##
proc addFoRCe {force msr T1 T2 T3 T4} {
    global atmcrd ts
    if {$T2!=$T4} {
        foreach {g1 g2 g3 g4} [dihedralgrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)] {}
        addforce $T1 [vecscale $g1 $force]
        addforce $T2 [vecscale $g2 $force]
        addforce $T3 [vecscale $g3 $force]
        addforce $T4 [vecscale $g4 $force] 
    } else {
        set Vecsub [vecsub $atmcrd($T1) $atmcrd($T2)]
        set force [expr $force / $msr]
        addforce $T1 [vecscale  $force $Vecsub]
        addforce $T2 [vecscale  [expr -$force] $Vecsub]
    }
} 
##------------------------------------------------------------##

proc addForceAngle {force T1 T2 T3} {
    global atmcrd ts
    foreach {g1 g2 g3} [anglegrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3)] {}
    addforce $T1 [vecscale $g1 $force]
    addforce $T2 [vecscale $g2 $force]
    addforce $T3 [vecscale $g3 $force]
}
# ############################################### 
# ##### Angle & Aihedral Force On Torsion ####### 
# ###############################################
proc AddForceB {const T1 T2} {
    global atmcrd
    set Vecsub [vecsub $atmcrd($T1) $atmcrd($T2)]
    addforce $T1 [vecscale  $const $Vecsub]
    addforce $T2 [vecscale  [expr -$const] $Vecsub]
}
##------------------------------------------------------------##

set pFile [open [format output_%03d/dat/Tor.%03d.dat $i_job $replica_id] w]
set qFile [open [format output_%03d/dat/img.%03d.dat $i_job $replica_id] w]
# set gFile [open [format output_%03d/dat/T2r.%03d.dat $i_job $replica_id] w]
set rFile [open [format output_%03d/dat/Prc.%03d.dat $i_job $replica_id] w]
fconfigure $pFile -buffering full -buffersize 1000000
# fconfigure $gFile -buffering full -buffersize 1000000
fconfigure $qFile -buffering line
fconfigure $rFile -buffering line

set PI 3.14159265358979
set deg2Rad  [expr $PI/180.0]
set deg2Rad2 [expr $deg2Rad*$deg2Rad / 2.0] ; # There is K/2 in harmonic potential
set oldTS -1
IndF $refID

##--------------------------------------##
source W117a.tcl
source F116a.tcl
source L359.tcl 
source L359.tcl 
source F116b.tcl 
source L366.tcl
source W117b.tcl 
source G358.tcl 
source F355.tcl 
# ######################################## 
# ############# Tcl Forces ############### 
# ######################################## 
proc calcforces {} {

    global kw pFile deg2Rad deg2Rad2 refList wr kw_twic Energy sec Run kwb kwb_twic Eng_b HelixList phiD   subsec  tiD i1 i2 i3 i4 j1 j2 j3 j4 k1 k2 k3 k4 l1 l2 l3 l4 m1 m2 m3 m4 Prc_362    kwB    Prc_HBond358   Prc_Bond BondToT   wrH rang_bond Kw_Bond  Prc_Ang Prc_Bnd phiBnd phiAng          aid_Ang Alpha     Prc_Bond360 Rc1 Rc2 Rc3 KwD BonD AngD Tor addFoRCe
    global atmcrd ts gFile rang IndexList cmForce  alpha2 rFile tFile table Prc_Helix Prc_Flip  wrB ranG aiD   Prc_113  PI Prc_Linear tor Linear tor step condF condA condT subsec0 tiA  alpha_UpBnd beta_UpBnd alpha_DnBnd beta_DnBnd KwBnd    phiD1 bond358  Kw_Bond360 Ang1 A366 HBond Prc_HBond alpha lDiff aid_Bnd Bond LinAngle addForceAngle   Rc4 Bd Ag  n1 n2 n3 n4
    global Prc_Flip force addForce alpha_Up    beta_Up    alpha_Dn    beta_Dn      oldTS IndF refID angleFixer DelPhi Eng qFile lable  Prc_PhiPsi Prc_Rc hFile addForceA RC phi362   alpha_Dn362    addFORCE stt aid_Ang221B aid_Bnd358   WrA   Ang366 AddForceB angleFixer2 Tor bond BoNd aid_Bnd   aid_Ang BnD Tor366 Tor6 Rc5 Rc6 Tag tiA358 tiA366 tiA359 omega1 omega2 omega3 omega4 omega5 omega6  omega7 omega8 omega9 omega10 omega11 omega12 omegaweight1 omegaweight2 omegaweight3 omegaweight4 omegaweight5 omegaweight6 omegaweight7 omegaweight8 omegaweight9 omegaweight10 omegaweight11 omegaweight12 beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10 beta11 beta12 OmegaWeight359 OmegaWeight366 OmegaWeight358 Alpha116a HBond116a Alpha116b HBond116b Alpha117a HBond117a Alpha117b HBond117b Alpha359 HBond359 Bond359 Ang359  Alpha366 HBond366 Bond366 Alpha358 HBond358 Alpha355 HBond355    aid_Ang366 aid_Ang359 aid_Ang358 US116a US359 US366 US358 US116b US117a US117b US366 US355 CV366 CV359 CV358 CV355 CV117b CV117a CV116a CV116b ColVar RCsub Exchange116a Exchange117a  Exchange117b Exchange355 Exchange358 Exchange359 Exchange366 Exchange116b aid_Ang117 OmegaWeight117 omega13 omega14 beta13 beta14 tiA117 Bound117a Bound117b Bound359 Bound358 Bound366 Bond358 AlphaRest EnergyBoundary OmegaWeight117b tiA117b omega15 omega16 omega17 beta15 beta16 beta17 Rcm

    
    set ts [getstep]
    loadcoords atmcrd
    

    
    if {$Tag == "R117a"} {
        US117a 
    } elseif {$Tag == "R116a"} {
        US116a
    } elseif {$Tag == "R359"}  {
        US359
    } elseif {$Tag == "R116b"} {
        US116b
    } elseif {$Tag == "R366"}  {
        US366 
    } elseif {$Tag == "R117b"} {
        US117b
    } elseif {$Tag == "R358"}  {
        US358
    } elseif {$Tag == "R355"}  {
        US355 
    }
    
return 
}

