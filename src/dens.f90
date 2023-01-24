
module modirtsre2

  implicit none
  integer,save ::ny,ns,nv,nea &
       ,maxmes,nobs &
       ,npmtot &
       ,nmescur &
       ,nvarxevt,nbevt,nrisqtot,nasso
  integer,save::nef_s,nvc_s,nef_B,nmu_v,nvc_v,nvc_u,nvc_err     
  double precision,dimension(:),allocatable,save::Y,minY,maxY
  double precision,dimension(:,:),allocatable,save ::X
  integer,dimension(:),allocatable,save::idg
  integer,dimension(:,:),allocatable,save::nmes
  !double precision,dimension(:),allocatable,save :: seqMC
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::b
  
end module modirtsre2




subroutine dens(re0,best0,Y0,X0,idg0 &
     ,ny0,ns0,nv0,nobs0,nmes0,npm0,nea0,zitr0 &
     ,ind0,nmescur0 &
     ,dens_res)

  use modirtsre2

  IMPLICIT NONE
  
  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0
  integer, intent(in)::ns0,nobs0,npm0,nea0
  !integer,intent(in)::idtrunc0,logspecif0,nbevt0
  !double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  !integer, dimension(ns0),intent(in)::Devt0
  !integer, dimension(nv0),intent(in)::idsurv0
  !integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  !double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
  double precision,dimension(2,ny0),intent(in)::zitr0
  integer, dimension(nv0),intent(in)::idg0
  integer,dimension(ns0,ny0),intent(in)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  double precision,dimension(nea0),intent(in)::re0
  double precision, dimension(npm0), intent(in) :: best0
  integer,intent(in)::ind0,nmescur0
  
  !Declaration des variables en sortie
  double precision,intent(out)::dens_res
  
  !Variables locales
  integer::jtemp,i,j,k,ktemp,yk,mi
  !integer::npmtot0
  double precision,external::pred
  
  !print*,"start"
  
  !print*,"best0",best0
  allocate(b(npm0))
  do k=1,npm0
     b(k)=best0(k)
  end do
  
  !print*,"#"
  !print*,"##"
  !print*,"### IN"
  !print*,"##"
  !print*,"#"
  
  !print*,"ELEMENTS EN ENTREE FORTRAN"
  !print*,"re0",re0
  !print*,"Y0",Y0
  !print*,"X0",X0
  !print*,"idg0",idg0
  !print*,"ny0",ny0
  !print*,"ns0",ns0
  !print*,"nv0",nv0
  !print*,"nobs0",nobs0
  !print*,"nmes0",nmes0
  !print*,"npm0",npm0
  !print*,"nea0",nea0
  !print*,"zitr0",zitr0
  !print*,"ind0",ind0
  !print*,"nmescur0",nmescur0
  !print*,"dens_res",dens_res
  !print*,"FIN ELEMENTS EN ENTREE FORTRAN"

  maxmes=0
  do i=1,ns0
     mi=sum(nmes0(i,:))
     if (mi.gt.maxmes) then
        maxmes=mi
     end if
  end do
  
  !rangeY
  allocate(minY(ny0),maxY(ny0))
  do k=1,ny0
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(2,k)
  end do
  
  allocate(Y(nobs0),X(nobs0,nv0) &
       ,idg(nv0),nmes(ns0,ny0))
       
  nbevt=0
  
  ny=ny0
  ns=ns0
  nv=nv0
  nobs=nobs0
  
  nmes=0
  Y=0.d0
  X=0.d0
  idg=0
  ktemp=0
  
  do k=1,nv
     idg(k)=idg0(k)
     
     jtemp=0
     do i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
              nmes(i,yk)=nmes0(i,yk)   !dim(nmes)=ns*ny    
              do j=1,nmes(i,yk)
                 jtemp=jtemp+1
                 Y(jtemp)=Y0(jtemp)
                 ktemp=ktemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           else
              do j=1,nmes(i,yk)
                 ktemp=ktemp+1
                 jtemp=jtemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           end if
        end do
     end do
  end do
  
  ! nb prms
  !npmtot0 = npm0
  
  
  ! creation des parametres
  
  nrisqtot=0
  
  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  
  nvarxevt = 0 !+ nvdepsurv

  nef_s=0
  do k=1,nv
     if (idg(k).eq.1) then
        nef_s = nef_s + 1
     end if
  end do
  
  nasso = nbevt*nea
  
  nvc_s = 1
  nef_B = ny
  nmu_v = ny
  nvc_v = ny
  nvc_u = ny 
  nvc_err = ny 
  
  nea = nvc_s + nvc_v + nvc_u ! nb RE
  
  npmtot = nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+nvc_u+nvc_err
  
  
  !print*,"deb pred"

  ! prediction des EAs
  dens_res = pred(re0,nea0,ind0,nmescur0)

  !print*,"fin pred"
  
  !print*,"dens_res",dens_res

  !1589 continue

  deallocate(Y,X,idg,nmes)

  deallocate(minY,maxY)
  
  deallocate(b)
  
  !print*,"#"
  !print*,"##"
  !print*,"### OUT"
  !print*,"##"
  !print*,"#"
  
  !print*,"end"
  !print*,"###"
  
  return
  
end subroutine dens



!----------------------------------------------!
!     INDIVIDUAL RANDOM EFFECTS PREDICTION     !
!----------------------------------------------!

double precision function pred(re,m,i,ncur)

  ! re = vecteur des effets aleatoires individuels
  ! m = nb re
  ! i = individu
  ! ncur = nb observations des individus precedents
  
  use modirtsre2

  implicit none

  integer::m,i,ncur
  double precision,dimension(m)::re
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer ::j,k,l
  integer ::sumMesYk,yk
  double precision,dimension(nef_s) :: b_s
  double precision,dimension(nef_s) :: X_s ! double precision,dimension(maxmes,nef_s) :: X_s
  double precision :: mu_s
  double precision,dimension(maxmes) :: Xtime
  double precision,dimension(maxmes) :: mu_v
  double precision,dimension(maxmes) :: rate
  double precision,dimension(maxmes) :: asympt_inf,asympt_sup
  double precision,dimension(maxmes,maxmes) :: var_Y, inv_var_Y
  double precision :: det_var_Y
  
  !double precision,dimension(maxmes) :: re_z
  !double precision,dimension(maxmes) :: re_v,re_u
  
  double precision,dimension(maxmes) :: shift, time_shifted
  double precision,dimension(maxmes) :: nu
  double precision,dimension(maxmes) :: denom,mu,u,mu_Y
  
  !double precision,dimension(nea,nea) ::Ut
  !double precision,dimension(npm) :: b
  double precision,dimension(npmtot)::b1
  !double precision,dimension(nxevt)::Xevt,bevt
  !double precision,dimension(maxval(nprisq))::brisq
  !double precision::basso 
  
  double precision ::Y4
  double precision,dimension(maxmes) :: Y1,Y2,Y3  !mu
  !double precision,dimension(nea)::ui,usim
  !double precision,dimension(nbevt)::risq,surv,surv0
  double precision::div,dens_Y  !,vrais_surv,varexpsurv
  !double precision::surv0_glob,surv_glob,fevt,easurv  
  !double precision::som_T0,som_Ti
  
  double precision,dimension(nea,nea) :: var_re, inv_var_re
  double precision :: det_var_re
  double precision ::Y4_re
  double precision,dimension(nea) :: mu_re,Y2_re,Y3_re
  double precision::div_re,dens_re
  
  pred=0.d0
  
  nmescur = ncur
  
  b1=0.d0
  do k=1,npmtot
     b1(k)=b(k)
  end do
  
  !print*, "b1", b1

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! creation de Y1 = vecteur des reponses observees
  Y1=0.d0
  sumMesYk = 0
  do yk=1,ny
      do j=1,nmes(i,yk)
          Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
      end do
      sumMesYk=sumMesYk+nmes(i,yk)
  end do !fin boucle yk
  
  !print*, "nmes(i,:)", nmes(i,:)
  !print*, "Y1", Y1
  
  !! Xtime
  Xtime=0.d0
  do j=1,sum(nmes(i,:))
      Xtime(j)=dble(X(nmescur+j,1))  
  end do
  
  !! rate
  rate=0.d0
  l=0
  do k=1,ny
     do j=1,nmes(i,k)
        l=l+1
        rate(l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+k)
     end do
  end do
  
  !! asymptotes
  asympt_inf=0.d0
  asympt_sup=0.d0
  l=0
  do k=1,ny
     do j=1,nmes(i,k)
        l=l+1
        asympt_inf(l)=minY(k)
        asympt_sup(l)=maxY(k)
     end do
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dens_Y=1.d0
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! vecteur des EAs
     
     !print*,"re",re
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie longitudinale ! f(Y|re)
     
     
     !!! time-shift
     
     !! s
     do j=1,sum(nmes(i,:))
        shift(j)=re(1)  
     end do
     !! t+s
     time_shifted = Xtime + shift
     
     !print*, "shift", shift
     !print*, "Xtime", Xtime
     !print*, "time_shifted", time_shifted
     
     !! v
     l=0
     do k=1,ny
        do j=1,nmes(i,k)
           l=l+1
           nu(l)=re(1+k)
        end do
     end do
     
     !print*, "nu", nu
     
     !! u
     l=0
     do k=1,ny
        do j=1,nmes(i,k)
           l=l+1
           u(l)=re(1+ny+k)
        end do
     end do
     
     !print*, "u", u
     
     
     !!! mu
     
     !! denominator
     denom = 0.d0
     do j=1,sum(nmes(i,:))
        denom(j) = 1 + exp( - rate(j)*time_shifted(j) )
        denom(j) = denom(j) ** nu(j)
     end do
     
     !print*, "rate", rate
     !print*, "denom", denom
     
     !! mu
     mu = 0.d0
     do j=1,sum(nmes(i,:))
        mu(j) = asympt_inf(j) + (asympt_sup(j) - asympt_inf(j)) / denom(j)
     end do
     
     !print*, "asympt_inf", asympt_inf
     !print*, "asympt_sup", asympt_sup
     !print*, "mu", mu
     
     !! mu_Y
     mu_Y = mu + u
     
     !print*, "mu_Y", mu_Y
     
    
     sumMesYk=0

     do yk =1,ny

          if(nmes(i,yk).gt.0) then
          
              !print*, "yk", yk
          
              !! var_Y = variance de Y|s,v,u
              ! avec inv_var_Y = matrice inverse de var_Y, comme var_Y diago, alors inverse = coefficient inverse
              ! avec det_var_Y = determinant de var_Y, comme var_Y diago, alors determinant = produit des coefficients diagonaux
              var_Y=0.d0
              inv_var_Y=0.d0
              det_var_Y=1
              do j=1,nmes(i,yk)
                  var_Y(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+nvc_u+yk)**2
                  !print*,"var_Y(j,j)",var_Y(j,j)
                  inv_var_Y(j,j) = 1/var_Y(j,j)
                  det_var_Y = det_var_Y * var_Y(j,j)
              end do
              
              ! calcul de la densite
              Y2=0.d0
              Y3=0.d0
              Y4=0.d0
              do j=1,nmes(i,yk)
                 Y2(j) = Y1(sumMesYk+j) - mu_Y(sumMesYk+j)
              end do
              
              !print*,"Y2",Y2
              
              Y3=matmul(inv_var_Y,Y2)
              
              !print*,"Y3",Y3
              
              Y4=DOT_PRODUCT(Y2,Y3)
              
              !print*,"Y4",Y4
              
              div = (dble(2*3.14159265)**(dble(nmes(i,yk))/2))*sqrt(det_var_Y) !exp(det_var_Y)
              
              !print*,"div",div
              
              !print*,"density",exp(-Y4/2.d0)/div
              
              dens_Y = dens_Y * exp(-Y4/2.d0)/div
              
          end if
        
          sumMesYk = sumMesYk+nmes(i,yk)
     end do ! fin boucle yk
     
     !print*,"dens_Y",dens_Y
  
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie survie deleted
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie effets aleatoires ! f(re)
  
     
     !! var_re = variance de s,v,u
     ! avec inv_var_re = matrice inverse de var_re, comme var_re diago, alors inverse = coefficient inverse
     ! avec det_var_re = determinant de var_re, comme var_re diago, alors determinant = produit des coefficients diagonaux
     var_re=0.d0
     inv_var_re=0.d0
     det_var_re=1
     do j=1,nea
       if(j.eq.1) var_re(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+1)**2  ! var_s
       if(j.gt.1) then
          if(j.le.(1+ny)) var_re(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+j-1)**2 ! var_v
          if(j.gt.(1+ny)) var_re(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+j-1-ny)**2 ! var_u
       end if
       !print*,"var_re(j,j)",var_re(j,j)
       inv_var_re(j,j) = 1/var_re(j,j)
       det_var_re = det_var_re * var_re(j,j)
     end do
     
     !print*,"var_re",var_re
     !print*,"inv_var_re",inv_var_re
     !print*,"det_var_re",det_var_re
     
     
     !! mu_re = moyenne de s,v,u
     mu_re=0.d0
     ! mu_s
     mu_s=0
     b_s=0.d0
     X_s=0.d0
     l=0
     do k=2,nv !colonne 1 = tps
        if (idg(k).ne.0) then
            l=l+1
            X_s(l)=dble(X(nmescur+1,k))  
            b_s(l)=b1(nrisqtot+nvarxevt+nasso+l)
         end if
     end do
     mu_s = dot_product(X_s,b_s) !matmul
     mu_re(1)=mu_s ! mu_s
     ! mu_v
     mu_v=0.d0
     l=0
     do k=1,ny
        l=l+1
        mu_v(l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+k)
     end do
     do j=1,nvc_v
       mu_re(1+j)=mu_v(j)
     end do
     ! mu_u
     !! moyennes = 0
     
     !print*,'mu_s',mu_s
     !print*,'mu_v',mu_v
     !print*,'mu_re',mu_re
     
     
     ! calcul de la densite
     dens_re=1.d0
     Y2_re=0.d0
     Y3_re=0.d0
     Y4_re=0.d0
     do j=1,nea
       Y2_re(j) = re(j) - mu_re(j)
     end do
     
     !print*,"Y2_re",Y2_re
              
     Y3_re=matmul(inv_var_re,Y2_re)
     
     !print*,"Y3_re",Y3_re
              
     Y4_re=DOT_PRODUCT(Y2_re,Y3_re)
     
     !print*,"Y4_re",Y4_re
              
     div_re = (dble(2*3.14159265)**(dble(nea)/2))*sqrt(det_var_re)
     
     !print*,"div_re",div_re
        
     dens_re = dens_re * exp(-Y4_re/2.d0)/div_re  
     !print*,"dens_re",dens_re
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  pred = dens_Y * dens_re !density
  pred = log(pred) !density in log
     
  !print*, "pred", pred
  
  if (pred.eq.-1.d9 .or. pred/pred.ne.1) then 
    pred = -1.d9
    !goto 541
  end if
  
  !654 continue

  return
  
end function pred  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!