       do ispin=1,Nspin         !Spin diagona
          do iorb=1,Norb        !Orbital diagonal part GF_0=(iw+mu)_aa-hloc_aa-Delta_aa
             do i=1,NL
                iw=xi*wm(i)
                impG0iw(ispin,ispin,iorb,iorb,i)= iw+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                impG0wr(ispin,ispin,iorb,iorb,i)= iw+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
          enddo
          do iorb=1,Norb         !Orbital non-diagonal part
             do jorb=iorb+1,Norb !GF_0=-hloc_ab-Delta_ab
                do i=1,NL
                   iw=xi*wm(i)
                   impG0iw(ispin,ispin,iorb,jorb,i)= -hloc(ispin,ispin,iorb,jorb)-delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                   impG0iw(ispin,ispin,jorb,iorb,i)= -hloc(ispin,ispin,jorb,iorb)-delta_bath(ispin,jorb,iorb,iw,dmft_bath)
                enddo
                do i=1,Nw
                   iw=cmplx(wr(i),eps)
                   impG0wr(ispin,ispin,iorb,jorb,i)= -hloc(ispin,ispin,iorb,jorb)-delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                   impG0wr(ispin,ispin,jorb,iorb,i)= -hloc(ispin,ispin,jorb,iorb)-delta_bath(ispin,jorb,iorb,iw,dmft_bath)
                enddo
             enddo
          enddo
       enddo
       !
       !                         !Get Sigma and G_0 by matrix inversions:
       do ispin=1,Nspin
          do i=1,NL
             invGimp = impGmats(ispin,ispin,:,:,i)
             impG0   = impG0iw(ispin,ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSmats(ispin,ispin,:,:,i) = impG0 - invGimp
             call matrix_inverse(impG0)
             impG0iw(ispin,ispin,:,:,i)=impG0
          enddo
          do i=1,Nw
             invGimp = impGreal(ispin,ispin,:,:,i)
             impG0   = impG0wr(ispin,ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSreal(ispin,ispin,:,:,i) = impG0 - invGimp
             call matrix_inverse(impG0)
             impG0wr(ispin,ispin,:,:,i)=impG0
          enddo
       enddo
       !
       !Print the impurity functions:
       do iorb=1,Norb
          do jorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
             call open_units(reg(suffix))
             do i=1,NL
                write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impGmats(ispin,ispin,iorb,jorb,i)),dreal(impGmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSmats(ispin,ispin,iorb,jorb,i)),dreal(impSmats(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impG0iw(ispin,ispin,iorb,jorb,i)),dreal(impG0iw(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Nw
                write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impGreal(ispin,ispin,iorb,jorb,i)),dreal(impGreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSreal(ispin,ispin,iorb,jorb,i)),dreal(impSreal(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impG0wr(ispin,ispin,iorb,jorb,i)),dreal(impG0wr(ispin,ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             call close_units()
          enddo
       enddo
       write(LOGfile,*)""
