       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,NL
                iw=xi*wm(i)
                fg0                        = iw+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSmats(ispin,ispin,iorb,iorb,i)= fg0 - one/impGmats(ispin,ispin,iorb,iorb,i)
                impG0iw(ispin,ispin,iorb,iorb,i) = one/fg0
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                fg0                        = wr(i)+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSreal(ispin,ispin,iorb,iorb,i)= fg0 - one/impGreal(ispin,ispin,iorb,iorb,i)
                impG0wr(ispin,ispin,iorb,iorb,i) = one/fg0
             enddo
          enddo
       enddo
       !
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))
          call open_units(reg(suffix))
          do i=1,NL
             write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impGmats(ispin,ispin,iorb,iorb,i)),dreal(impGmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impSmats(ispin,ispin,iorb,iorb,i)),dreal(impSmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impG0iw(ispin,ispin,iorb,iorb,i)),dreal(impG0iw(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          do i=1,Nw
             write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impGreal(ispin,ispin,iorb,iorb,i)),dreal(impGreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impSreal(ispin,ispin,iorb,iorb,i)),dreal(impSreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impG0wr(ispin,ispin,iorb,iorb,i)),dreal(impG0wr(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          call close_units
       enddo
