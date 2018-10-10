c***********************************************************
c display format  1                                        *
c implicit scheme 2                                        *
c 0.002
c******************DeltT~1s,10s*****************************
c****************************************Zhu 18/9/2016******
      Module Tsd
      Implicit none
c-----Module shape_4s8n----------------
c-------nodes,dispacement--------------
c-------coordanary---------------------
      Integer(kind=4) :: imax,jmax
c----节点总数，单元数，单元节点，单节点位移数，位移边界个数
c----承受外载数，物性参数个数，多点约束，总位移数
      Integer(kind=4) :: nn,ne,nen,ndn,nd,nl,npr,nmpc,nq,nbw
      Integer(kind=4) :: NOC(2000,8),NU(800)
      Real(kind=8) :: rad
      Real(kind=8) :: X(4000,2),DT(4000),DT0(2000),PM(6),D(4,4),DD(4,4),
     1    B(4,16),DB(4,16),SE(16,16),S(8000,800),TL(16),DPL(16),DTL(16),
     1    Q(16),U(300),F(8000),
c     1    X0(2000),Y0(2000),T0(2000),
     1    T01(4000),T00(4000),FLAG(4000),DDT(4000)
      Real(kind=8),Allocatable::X0(:),Y0(:),T0(:)
c----set T00(:) aimming to conserve the initial te_data to compute dte
c-----Module Dis_den-------------------
c-----文献中的gshear对应于本程序中的mater中的弹性参数
      Real(kind=8) :: bv0,strainhard,acten,strexp,matecon,boltcon,kv0
c-----Module tymath--------------------      
c------implicit scheme 2---------------
      Real(kind=8) :: j2,efs,g,f1(4),f2,strdev(4),dgdj2,dgdnm
      Real(kind=8) :: dj2dsigma(4),df1dsigma(4,4),df1dnm(4),df1dt(4),
     1    df2dsigma(4),df2dnm,df2dt,hn(4,4),talpha(4)
c-----talpha(:) for considering the anisotropy of therm-deformation
c-----Numerical integration------------
      Integer(kind=4) :: NIK
      Real(kind=8) :: DJ,TJ11,TJ12,TJ21,TJ22
      Real(kind=8) :: NIX(4,2)
c-----Numerical scheme----------------- 
c-----Module tsd,时间，隐式因子，绝对温度---------------------
      Real(kind=8) :: Tpass,Ttol,Theta,DeltT,tem,Tpass2,Ttol2
      Integer(kind=4)::NT,NT0,NT2
c-----Module transform,函数功能转换,功能控制
      Integer(kind=4)::issol
c-----那么TDENS,TSTRESS里面装的是n时刻到n+1时刻过程中位错以及应力的增量
c-----(单元，积分点)格式
      Real(kind=8) :: DENS(4000),TDENS(4000),STRESS(4000,4),
     1    TSTRESS(4000,4),STR(4000,4,4),VMSTRESS(4000),INDENS(4000,4),
     1    TSTR(4000,4,4),TINDENS(4000,4)
      Contains
      Function thermalpha(local)
      Implicit None
      Real(kind=8)::thermalpha,local
      Real(kind=8)::a1=-5.54e-6,a2=-1.39e-4,a3=-1.39e-4
      Real(kind=8)::t1=3358.8,t2=42.3,t3=42.3
      Real(kind=8)::alpha0=8.83e-6
      thermalpha=a1*exp(-local/t1)+a2*exp(-local/t2)+
     1           a3*exp(-local/t3)+alpha0
      Return
      End Function thermalpha

      End Module Tsd

      Program tsd_1
      Use Tsd
      Integer(kind=4)::i,j,cnum
      integer(kind=4) :: state
      c=0.8660254
c      Ttol=1*3600
c      Ttol=4*1800
      Ttol=12*1800
      Tpass=0
      Tpass2=0
      NT=0
c      NT0=1
      NT0=0
      NT2=0
c      Theta=1.0
      Theta=0.0
      DeltT=0.002
      write(*,'(3f12.4)')Ttol,DeltT,Theta
      cnum=int(Ttol/DeltT)+1
      write(*,*)'the total steps number:',cnum
      state=0
      write(*,*)'whether there is input data',state
      call initmod()
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      If(state.eq.1)then
      Call readinput()
      write(*,99)NT0,Tpass,Tpass2,DeltT
      else
      write(*,100)NT,Tpass
      Call initsim()
      Endif
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'(6es13.5)')VMSTRESS(49),DENS(49),(STRESS(49,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(565),DENS(565),
     1(STRESS(565,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1009),DENS(1009),
     1(STRESS(1009,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1727),DENS(1727),
     1(STRESS(1727,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1751),DENS(1751),
     1(STRESS(1751,i),i=1,4)
c      do NT=1,int(Ttol/DeltT)+1
c      Do NT=1,cnum
      If(Tpass2.eq.0)then
      Do NT=NT0+1,cnum
      If(NT.gt.10000)DeltT=0.01
c      If(NT.gt.100000)DeltT=0.01
c      Do NT=1,20
      Tpass=Tpass+DeltT
      write(*,100)NT,Tpass
99    format(1x,'time num NT0=',i8,' the initial time Tinitial=',f12.4)
100   format(1x,'time num NT=',i8,' tpass time =',f12.3)
      call solstrt() 
      write(*,'(6es13.5)')VMSTRESS(49),DENS(49),(STRESS(49,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(565),DENS(565),
     1(STRESS(565,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1009),DENS(1009),
     1(STRESS(1009,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1727),DENS(1727),
     1(STRESS(1727,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1751),DENS(1751),
     1(STRESS(1751,i),i=1,4)
 
      If(Tpass.gt.Ttol)exit
      enddo 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NT=cnum
c      write(*,*)'initial program'
c      write(*,*)'NT the second stage'
      else
      do NT2=1,10000
      Tpass2=Tpass2+DeltT
      write(*,101)NT2,Tpass2
101   format(1x,'residual time num NT2=',i8,'total time Tpass2=',f12.4) 
      call solstrt()
      write(*,'(6es13.5)')VMSTRESS(49),DENS(49),(STRESS(49,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(565),DENS(565),(STRESS(565,i),
     1 i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1009),DENS(1009),(STRESS(1009,i),
     1 i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1727),DENS(1727),
     1(STRESS(1727,i),i=1,4)
      write(*,'(6es13.5)')VMSTRESS(1751),DENS(1751),
     1(STRESS(1751,i),i=1,4)
      enddo
      Endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      stop
      End Program

      Subroutine readinput()
      Use Tsd
      Implicit none
      Integer(kind=4)::i,j,k,num
      open(10,file='input')
      read(10,'(3f12.4,i8)')Ttol,Tpass,Tpass2,DeltT,Theta,NT0
      Read(10,*)
      Read(10,*)
      do j=1,jmax
      do i=1,imax
      num=(j-1)*(3*imax-1)+2*i-1
      read(10,001)X(num,1),X(num,2),T00(num),T01(num),VMSTRESS(num),
     1DENS(num),(STRESS(num,k),k=1,4)
      enddo
      enddo
      read(10,*)
      do i=1,ne
      do j=1,4
      read(10,002)(STR(i,j,k),k=1,4),INDENS(i,j)
      enddo
      enddo
      close(10)
001   format(2f10.5,8es15.5)
002   format(5es15.5)
      Write(*,*)'Call readinput'
      return
      End Subroutine

      Subroutine chtem()
      Use Tsd
      Implicit none
      Integer(kind=4)::i,j,num,nod
      Real(kind=8)::Tm,Taver
      Real(kind=8),Allocatable::DTt(:)
      Real(kind=8)::T_MAX=0,T_MIN=1e9
      Allocate(DTt(size(DT0)))
      Tm=800
      Do i=1,imax
      Do j=1,jmax
      num=(j-1)*(3*imax-1)+2*i-1
      nod=i+(j-1)*imax
c      T01(num)=(T0(nod)-Tm)/Ttol*(Ttol-NT*DeltT)+Tm
      T01(num)=(T0(nod)-Tm)/Ttol*(Ttol-Tpass)+Tm
c      T01(num)=T00(num)
      If(T01(num).lt.T_MIN)T_MIN=T01(num)
      If(T01(num).gt.T_MAX)T_MAX=T01(num)
      Enddo
      Enddo
      Taver=(T_MAX+T_MIN)/2.0
      Do i=1,imax
      Do j=1,jmax
      num=(j-1)*(3*imax-1)+2*i-1
      nod=i+(j-1)*imax
      DTt(nod)=T01(num)-Taver
      Enddo
      Enddo
      Do i=1,imax
      Do j=1,jmax
      num=(j-1)*(3*imax-1)+2*i-1
      nod=i+(j-1)*imax
      DDT(num)=DTt(nod)-DT(num)
      Enddo
      Enddo
      Do i=1,imax
      Do j=1,jmax
      num=(j-1)*(3*imax-1)+2*i-1
      nod=i+(j-1)*imax
      DT(num)=DTt(nod)
      Enddo
      Enddo
      Return
      End subroutine

      Subroutine initmod
c----read in temperature,code and record model information
      Use Tsd
      Implicit none
      Integer(kind=4)::i,j,inum,jnum,nno,sign,tmp,num,nod
      Integer(kind=4)::sum=0
      Character*60 dummy
      Real(kind=8)::T_MIN=1E9,T_MAX=0,T_AVER,Y_MAX
      Real(kind=8)::nmin,nmax,ntmp
      Namelist /paraden/ bv0,strainhard,acten,strexp,matecon,boltcon,kv0
      Real(kind=8)::c=0.57735026919
      Open(10,file='nearcrys.dat')
c1000	format(1x,'zone t= "crystal",i=',i3,' , j=',i3,' , f= point')
      read(10,*)
      read(10,*)
      read(10,*)
      Read(10,*)inum,jnum
      nno=inum*jnum
      Allocate(X0(nno),Y0(nno),T0(nno))
      
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      Do i=1,nno
c      Read(10,*)sign,X0(i),Y0(i),T0(i)
c      Read(10,*)X0(i),Y0(i),T0(i)
c      Enddo
      Do j=1,jnum
      Do i=1,inum
      num=(jnum-j)*inum+i
      Read(10,*)X0(num),Y0(num),T0(num)
      Enddo
      Enddo
      Close(10)
      Open(11,file='mater.dat')
      Read(11,'(a)')dummy
      Read(11,*) (PM(i),i=1,6)
      Close(11)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      T0=T0-200
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Y_MAX=Y0(1)
      Do i=1,nno
      X0(i)=X0(i)-X0(1)
      Y0(i)=abs(Y0(i)-Y_MAX)
      If(T0(i).lt.T_MIN) T_MIN=T0(i)
      If(T0(i).gt.T_MAX) T_MAX=T0(i)
      Enddo
      T_AVER=(T_MAX+T_MIN)/2
      Do i=1,nno
      DT0(i)=T0(i)-T_AVER
      Enddo
c-----DT0原始温度不均匀度,T0原始绝对温度(没有全局编号时的原始参数)-------
c---------the total message--
      imax=inum
      jmax=jnum
c------nn:total nodes,ne:total elements
c----节点总数，单元数，单元节点，单节点位移数，位移边界个数
c----承受外载数，物性参数个数，多点约束，总位移数
      nn=jmax*(2*imax-1)+imax*(jmax-1)
      ne=(imax-1)*(jmax-1)
      nen=8
      ndn=2
      nd=(2*imax-1)*2+2*jmax-2
      nl=0
      npr=7
      nmpc=0
      nq=ndn*nn
c--------nodes message-------
      Do j=1,jmax-1
      Do i=1,imax-1
      tmp=(imax-1)*(j-1)+i
      NOC(tmp,1)=(j-1)*(3*imax-1)+2*i-1
      NOC(tmp,2)=(j-1)*(3*imax-1)+2*i+1
      NOC(tmp,3)=NOC(tmp,2)+3*imax-1
      NOC(tmp,4)=NOC(tmp,1)+3*imax-1
      Enddo
      Enddo

      Do j=1,jmax
      Do i=1,imax
      num=(j-1)*(3*imax-1)+2*i-1
      nod=i+(j-1)*imax
      X(num,1)=X0(nod)
      X(num,2)=Y0(nod)
      DT(num)=DT0(nod)
      T01(num)=T0(nod)
      T00(num)=T0(nod)
c----DT,T01分别是经过全局标号后的每个节点上的温度不均匀度及绝对温度
      Enddo
      Enddo

      Do j=1,jmax-1
      Do i=1,imax-1
      tmp=(imax-1)*(j-1)+i
      NOC(tmp,5)=NOC(tmp,1)+1
      NOC(tmp,7)=NOC(tmp,4)+1
      NOC(tmp,8)=NOC(tmp,1)+(imax-i)*2+i
      NOC(tmp,6)=NOC(tmp,8)+1
      X(NOC(tmp,5),1)=0.5*(X(NOC(tmp,1),1)+X(NOC(tmp,2),1))
      X(NOC(tmp,5),2)=X(NOC(tmp,1),2)
      X(NOC(tmp,7),1)=0.5*(X(NOC(tmp,4),1)+X(NOC(tmp,3),1))
      X(NOC(tmp,7),2)=X(NOC(tmp,4),2)
      X(NOC(tmp,8),1)=X(NOC(tmp,1),1)
      X(NOC(tmp,8),2)=0.5*(X(NOC(tmp,1),2)+X(NOC(tmp,4),2))
      X(NOC(tmp,6),1)=X(NOC(tmp,2),1)
      X(NOC(tmp,6),2)=0.5*(X(NOC(tmp,2),2)+X(NOC(tmp,3),2))
      Enddo
      Enddo
c-----sum 记录位移条件总数，NU记录边界条件节点位移顺序(位移向量中的位置)
c-----U(:)记录边界条件节点的位移
c-----边界条件是顶端固定及中间对称轴条件
c***************************************
      do i=1,imax-2
      NU(sum+1)=2*NOC(i,1)-1
      NU(sum+2)=2*NOC(i,1)
      NU(sum+3)=2*NOC(i,5)-1
      NU(sum+4)=2*NOC(i,5)
      sum=sum+4
      enddo
c      do i=1,imax-2
c      NU(sum+1)=2*NOC(i,1)
c      NU(sum+2)=2*NOC(i,5)
c      sum=sum+2
c      enddo
c      NU(sum+1)=2*NOC(imax-1,1)
c      NU(sum+1)=2*NOC(imax-1,5)
c      NU(sum+3)=2*NOC(imax-1,2)
c      sum=sum+3

      NU(sum+1)=2*NOC(imax-1,1)-1
      NU(sum+2)=2*NOC(imax-1,1)
      NU(sum+3)=2*NOC(imax-1,5)-1
      NU(sum+4)=2*NOC(imax-1,5)
      NU(sum+5)=2*NOC(imax-1,2)-1
      NU(sum+6)=2*NOC(imax-1,2)
      sum=sum+6

      do i=1,(jmax-2)*(imax-1)+1,imax-1
      NU(sum+1)=2*NOC(i,8)-1
      NU(sum+2)=2*NOC(i,4)-1
      sum=sum+2
      enddo
c      do i=imax-1,(jmax-1)*(imax-1),imax-1
c      NU(sum+1)=2*NOC(i,6)-1
c      NU(sum+2)=2*NOC(i,6)
c      NU(sum+3)=2*NOC(i,3)-1
c      NU(sum+4)=2*NOC(i,3)
c      sum=sum+4
c      enddo

      do i=1,nd
      U(i)=0
      enddo
c----计算该带宽----------
      do i=1,ne
      nmin=NOC(i,1)
      nmax=NOC(i,1)
      do j=2,4
      if(nmin.gt.NOC(i,j))nmin=NOC(i,j)
      if(nmax.lt.NOC(i,j))nmax=NOC(i,j)
      enddo
      ntmp=ndn*(nmax-nmin+1)
      if(nbw.lt.ntmp)nbw=ntmp
      enddo
c---数值积分点，注意相对位置
      NIX(1,1)=-c
      NIX(1,2)=-c
      NIX(2,1)=c
      NIX(2,2)=-c
      NIX(3,1)=c
      NIX(3,2)=c
      NIX(4,1)=-c
      NIX(4,2)=c
c---subroutine parain()
      Open(3,file='paran.dat',status='old')
      Read(3,nml=paraden)
      Close(3)
      Open(10,file='element.dat')
      Do i=1,ne
      Write(10,'(9i10)')i,(NOC(i,j),j=1,8)
      Enddo
      Close(10)
      Open(10,file='nodes.dat')
      Do i=1,nn
      Write(10,'(i6,2f10.5)')i,X(i,1),X(i,2)
      Enddo
      Close(10)
      open(10,file='nu.dat')
      do i=1,nd
      write(10,'(2i6,f10.1)')i,NU(i),U(i)
      enddo
      close(10)
      write(*,*)'bandwidth complete,nbw=',nbw
      write(*,*)'init function works!'
      Return
      End Subroutine

      Subroutine initsim()
      Use Tsd
      Implicit None
      Integer(kind=4) :: i,j,k,eind
      Integer(kind=4) :: iin,ii
      Real(kind=8)::c,c1,c2
      c=0.8660254
      issol=0
c---element load TL(:) and physical model total load F(:) init
      Do i=1,nq
      F(i)=0
      Do j=1,nbw
      S(i,j)=0
      Enddo
      Enddo

      Do eind=1,ne
      Do i=1,16
      Do j=1,16
      SE(i,j)=0
      Enddo
      TL(i)=0
      Enddo
      Do NIK=1,4
      Call dbmat_0(eind)
      Call elek()
      Call temstr(eind)
      Enddo
      Call eletoall(eind) 
      Enddo
      Call bansol()
c----all matrixes init first
      Do i=1,nn
      FLAG(i)=0
      VMSTRESS(i)=0
      DENS(i)=2.0e+6
      TDENS(i)=0
      Do j=1,4
      STRESS(i,j)=0
      TSTRESS(i,j)=0
      Enddo
      Enddo

      Do i=1,ne
      Do j=1,4
      INDENS(i,j)=2.0e+6
      TINDENS(i,j)=0
      Do k=1,4
      STR(i,j,k)=0
      TSTR(i,j,k)=0
      Enddo
      Enddo
      Enddo
c---computing the stress with the known displace
      issol=1
      Do eind=1,ne
      Do NIK=1,4
      Call dbmat_0(eind)
      Enddo
      Enddo
c---计算插值节点应力
      Do eind=1,ne
      Do i=1,4    
      FLAG(NOC(eind,i))=FLAG(NOC(eind,i))+1
      Enddo
      Do i=1,4
      STRESS(NOC(eind,1),i)=STRESS(NOC(eind,1),i)+(1+c)*STR(eind,1,i)+
     1  (1-c)*STR(eind,3,i)-0.5*(STR(eind,2,i)+STR(eind,4,i))
      STRESS(NOC(eind,2),i)=STRESS(NOC(eind,2),i)+(1+c)*STR(eind,2,i)+
     1  (1-c)*STR(eind,4,i)-0.5*(STR(eind,1,i)+STR(eind,3,i))
      STRESS(NOC(eind,3),i)=STRESS(NOC(eind,3),i)+(1-c)*STR(eind,1,i)+
     1  (1+c)*STR(eind,3,i)-0.5*(STR(eind,2,i)+STR(eind,4,i))
      STRESS(NOC(eind,4),i)=STRESS(NOC(eind,4),i)+(1-c)*STR(eind,2,i)+
     1  (1+c)*STR(eind,4,i)-0.5*(STR(eind,1,i)+STR(eind,3,i))
      Enddo
      Enddo
c---对同一个节点进行平均，根据的是FLAG(NOC(eind,i)),表明同一个节点计算了几次
      Do j=1,jmax
      Do i=1,imax
      k=(j-1)*(3*imax-1)+2*i-1
      Do ii=1,4
      STRESS(k,ii)=STRESS(k,ii)/FLAG(k)
      Enddo
      c1=STRESS(k,1)+STRESS(k,2)+STRESS(k,4)
      c2=STRESS(k,1)*STRESS(k,2)+STRESS(k,2)*STRESS(k,4)+
     1STRESS(k,4)*STRESS(k,1)-STRESS(k,3)*STRESS(k,3)
      VMSTRESS(k)=sqrt(c1*c1-3*c2) 
      Enddo
      Enddo
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      if(NT2.ne.0)then
      Call zppostback()
c      endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      End Subroutine

      Subroutine bansol()
c---半带宽方程高斯求解过程
c---计算出来的是节点位移
      use Tsd
      Implicit None
      Integer(kind=4)::i,j,k,ii,jj,nk
      Integer(kind=4)::ln,uln,eln
      Real(kind=8)::c
c      Allocatable, Real(kind=8)::S1(:,:),F1(:)
c	character(len=50)::nfile1
c---'set 1 method' ensure the special nodes'displacment 0
c record the initl stiffness matrix
c Only when the F(:) don't vary and reach the staady, can the next step
c simulating begin.
c The simple method which doesn't bring many changes is to restore the
c S(:,:).
c It's true that first we need to change S(:,:) through DENS(:,:) then
c recompulate S(:,:) and F(:,:),DENS(:,:).
c      Allocate(S1(size(S,1),size(S,2)),F1(size(F,1)))
c      do i=1,size(S,1)
c      do j=1,size(S,2)
c      S1(i,j)=S(i,j)
c      enddo
c      enddo
      Do i=1,nd
      ln=NU(i)
      uln=ln-1
      If(uln.lt.1)uln=1
      eln=ln-nbw
      If(eln.lt.1)eln=1
      Do j=uln,eln,-1
      S(j,ln-j+1)=0
      Enddo
      Do j=2,nbw
      S(ln,j)=0
      Enddo
c---给定位移为0
      S(ln,1)=1
      F(ln)=0
      Enddo
c---forward elimination
      Do k=1,nq-1
      nk=nq-k+1
      If(nk.gt.nbw)nk=nbw
      Do i=2,nk
      c=S(k,i)/S(k,1)
      ii=k+i-1
      Do j=i,nk
      jj=j-i+1
      S(ii,jj)=S(ii,jj)-c*S(k,j)
      Enddo
      F(ii)=F(ii)-c*F(k)
      Enddo
      Enddo
c---back substitution
      F(nq)=F(nq)/S(nq,1)
      Do ii=1,nq-1
      i=nq-ii
      c=1/S(i,1)
c---首项系数为1
      F(i)=c*F(i)
      nk=nq-i+1
      If(nk.gt.nbw)nk=nbw
      Do j=2,nk
      F(i)=F(i)-c*S(i,j)*F(i+j-1)
      Enddo
      Enddo
      Return
      End Subroutine

      Subroutine solstrt()
      Use Tsd
      Implicit None
      Integer(kind=4)::eind
      Integer(kind=4)::i,j,k,ii
      Real(kind=8)::c,c1,c2
      c=0.8660254
      issol=0
c----gain improvement of variables
c----the right of equations init
      Call chtem()
      Do i=1,nq
      F(i)=0
      Do j=1,nbw
      S(i,j)=0
      Enddo
      Enddo
      Do eind=1,ne
      Do i=1,16
      Do j=1,16
      SE(i,j)=0
      Enddo
      TL(i)=0
      Enddo
      Do NIK=1,4
      Call dbtmat_2(eind)
      Call elek()
      Call temstr_d(eind)
      Enddo
      Call eletoall(eind)
      Enddo
      Call bansol()
      Do i=1,nn
      FLAG(i)=0
      TDENS(i)=0
      Do j=1,4
      TSTRESS(i,j)=0
      Enddo
      Enddo

      Do i=1,ne 
      Do j=1,4
      TINDENS(i,j)=0
      Do k=1,4
      TSTR(i,j,k)=0
      Enddo
      Enddo
      Enddo

      issol=1
      Do eind=1,ne
      Do NIK=1,4
      Call dbtmat_2(eind)
      Enddo
      Enddo

      Do eind=1,ne
      Do i=1,4
      FLAG(NOC(eind,i))=FLAG(NOC(eind,i))+1
      Enddo
      Do i=1,4
      TSTRESS(NOC(eind,1),i)=TSTRESS(NOC(eind,1),i)+(1+c)*TSTR(eind,1,i)
     1+(1-c)*TSTR(eind,3,i)-0.5*(TSTR(eind,2,i)+TSTR(eind,4,i))
      TSTRESS(NOC(eind,2),i)=TSTRESS(NOC(eind,2),i)+(1+c)*TSTR(eind,2,i)
     1+(1-c)*TSTR(eind,4,i)-0.5*(TSTR(eind,1,i)+TSTR(eind,3,i))
      TSTRESS(NOC(eind,3),i)=TSTRESS(NOC(eind,3),i)+(1-c)*TSTR(eind,1,i)
     1+(1+c)*TSTR(eind,3,i)-0.5*(TSTR(eind,2,i)+TSTR(eind,4,i))
      TSTRESS(NOC(eind,4),i)=TSTRESS(NOC(eind,4),i)+(1-c)*TSTR(eind,2,i)
     1+(1+c)*TSTR(eind,4,i)-0.5*(TSTR(eind,1,i)+TSTR(eind,3,i))
      Enddo
      TDENS(NOC(eind,1))=TDENS(NOC(eind,1))+(1+c)*TINDENS(eind,1)
     1+(1-c)*TINDENS(eind,3)-0.5*(TINDENS(eind,2)+TINDENS(eind,4))
      TDENS(NOC(eind,2))=TDENS(NOC(eind,2))+(1+c)*TINDENS(eind,2)
     1+(1-c)*TINDENS(eind,4)-0.5*(TINDENS(eind,1)+TINDENS(eind,3))
      TDENS(NOC(eind,3))=TDENS(NOC(eind,3))+(1-c)*TINDENS(eind,1)
     1+(1+c)*TINDENS(eind,3)-0.5*(TINDENS(eind,2)+TINDENS(eind,4))
      TDENS(NOC(eind,4))=TDENS(NOC(eind,4))+(1-c)*TINDENS(eind,2)
     1+(1+c)*TINDENS(eind,4)-0.5*(TINDENS(eind,1)+TINDENS(eind,3))
      Enddo

      Do j=1,jmax
      Do i=1,imax
      k=(j-1)*(3*imax-1)+2*i-1
      Do ii=1,4
      TSTRESS(k,ii)=TSTRESS(k,ii)/FLAG(k)
      STRESS(k,ii)=STRESS(k,ii)+TSTRESS(k,ii)
      Enddo
      TDENS(k)=TDENS(k)/FLAG(k)
      DENS(k)=DENS(k)+TDENS(k)
      c1=STRESS(k,1)+STRESS(k,2)+STRESS(k,4)
      c2=STRESS(k,1)*STRESS(k,2)+STRESS(k,2)*STRESS(k,4)+
     1STRESS(k,4)*STRESS(k,1)-STRESS(k,3)*STRESS(k,3)
      VMSTRESS(k)=sqrt(c1*c1-3*c2) 
      Enddo
      Enddo
c      if(NT2.ne.0)then
      Call zppostback()
c      endif
      End Subroutine
c D matrix derive from ZZB content
c ~D matrix derive from D, vary with the time
c db matrix is S matrix
c dmat0----->>>>dmatrix
      Subroutine dmat_0()
c----原始的刚度矩阵，只是引入位错模型后刚度矩阵发生变化
      Use Tsd
      Implicit none
      Integer(kind=4)::i,j
      Real(kind=8)::c11,c12,c13,c33,c44
      c11=PM(1)
      c12=PM(2)
      c13=PM(3)
      c33=PM(4)
      c44=PM(5)
      Do i=1,4
      Do j=1,4
      D(i,j)=0
      Enddo
      Enddo
      D(1,1)=c11
      D(1,2)=c13
      D(1,4)=c12
      D(2,1)=D(1,2)
      D(2,2)=c33
      D(2,4)=c13
      D(3,3)=c44
      D(4,1)=D(1,4)
      D(4,2)=D(2,4)
      D(4,4)=c11
c -----H matrix----------------------        
      Return
      End subroutine

      Subroutine dmat_01(eind)
c----consider the temperature's influence
      Use Tsd
      Implicit none
      Integer(kind=4)::i,j,eind
      Real(kind=8)::c11,c12,c13,c33,c44
      Real(kind=8)::T_in,T_point(4)
      T_in=0
      Do i=1,4
      T_point(i)=T01(NOC(eind,i))
      T_in=T_in+T_point(i)
      Enddo
      T_in=T_in/4.0
      c11=PM(1)*exp(-9.4e-5*(T_in-298.15))
      c12=PM(2)*exp(-9.4e-5*(T_in-298.15))
      c13=PM(3)*exp(-9.4e-5*(T_in-298.15))
      c33=PM(4)*exp(-9.4e-5*(T_in-298.15))
      c44=PM(5)*exp(-9.4e-5*(T_in-298.15))
      Do i=1,4
      Do j=1,4
      D(i,j)=0
      Enddo
      Enddo
      D(1,1)=c11
      D(1,2)=c13
      D(1,4)=c12
      D(2,1)=D(1,2)
      D(2,2)=c33
      D(2,4)=c13
      D(3,3)=c44
      D(4,1)=D(1,4)
      D(4,2)=D(2,4)
      D(4,4)=c11
c -----transient DD matrix-----------
c -----H matrix---------------------- 
      Return
      End subroutine

      Subroutine dbmat_0(eind)
c---原始的应力函数矩阵
c---此函数是有功能转换的issol
      Use Tsd
      Implicit None
      Integer(kind=4)::eind
      Integer(kind=4):: i,j,k,iin,ii
      Integer(kind=4):: NODES(8)
      Real(kind=8)::dte,ddte
      Real(kind=8)::pi,xi,eta,c,c1,alpha
      Real(kind=8)::SH(8),SHX(8),SHY(8),SHT(4),R(8),Z(8)
      Real(kind=8)::AMAT(3,4),GMAT(4,16)
      pi=3.1415926
      xi=NIX(NIK,1)
      eta=NIX(NIK,2)
      Do i=1,8
      NODES(i)=NOC(eind,i)
      Enddo
      Do i=1,8
      R(i)=X(NODES(i),1)
      Z(i)=X(NODES(i),2)
      Enddo
c---shape function-----
      SH(1)=-0.25*(1-xi)*(1-eta)*(1+xi+eta)
      SH(2)=-0.25*(1+xi)*(1-eta)*(1-xi+eta)
      SH(3)=-0.25*(1+xi)*(1+eta)*(1-xi-eta)
      SH(4)=-0.25*(1-xi)*(1+eta)*(1+xi-eta)
      SH(5)=0.5*(1-xi*xi)*(1-eta)
      SH(6)=0.5*(1+xi)*(1-eta*eta)
      SH(7)=0.5*(1-xi*xi)*(1+eta)
      SH(8)=0.5*(1-xi)*(1-eta*eta)
c---shape derivative function
      SHX(1)=0.25*(1-eta)*(2*xi+eta)
      SHY(1)=0.25*(1-xi)*(2*eta+xi)
      SHX(2)=0.25*(1-eta)*(2*xi-eta)
      SHY(2)=0.25*(1+xi)*(2*eta-xi)
      SHX(3)=0.25*(1+eta)*(2*xi+eta)
      SHY(3)=0.25*(1+xi)*(2*eta+xi)
      SHX(4)=0.25*(1+eta)*(2*xi-eta)
      SHY(4)=0.25*(1-xi)*(2*eta-xi)
      SHX(5)=-xi*(1-eta)
      SHY(5)=-0.5*(1-xi*xi)
      SHX(6)=0.5*(1-eta*eta)
      SHY(6)=-eta*(1+xi)
      SHX(7)=-xi*(1+eta)
      SHY(7)=0.5*(1-xi*xi)
      SHX(8)=-0.5*(1-eta*eta)
      SHY(8)=-eta*(1-xi)
c---formation of jacobian
      TJ11=0
      TJ12=0
      TJ21=0
      TJ22=0
      Do i=1,8
      TJ11=TJ11+SHX(i)*R(i)
      TJ12=TJ12+SHX(i)*Z(i)
      TJ21=TJ21+SHY(i)*R(i)
      TJ22=TJ22+SHY(i)*Z(i)
      Enddo
c---determinant of jacobian
      DJ=TJ11*TJ22-TJ12*TJ21
c---A,G matrix ,geometry function
      Do i=1,3
      Do j=1,4
      AMAT(i,j)=0
      Enddo
      Enddo
      AMAT(1,1)=TJ22/DJ
      AMAT(1,2)=-TJ12/DJ
      AMAT(2,3)=-TJ21/DJ
      AMAT(2,4)=TJ11/DJ
      AMAT(3,1)=-TJ21/DJ
      AMAT(3,2)=TJ11/DJ
      AMAT(3,3)=TJ22/DJ
      AMAT(3,4)=-TJ12/DJ
      Do i=1,4
      Do j=1,16
      GMAT(i,j)=0
      Enddo
      Enddo
      Do i=1,8
      GMAT(1,2*i-1)=SHX(i)
      GMAT(2,2*i-1)=SHY(i)
      GMAT(3,2*i)=SHX(i)
      GMAT(4,2*i)=SHY(i)
      Enddo
c---任意一点处实际温度tem以及半径rad的表示
      rad=0
      Do i=1,8
      rad=rad+SH(i)*R(i)
      enddo
c---B matrix------------------
      Do i=1,3
      Do j=1,16
      B(i,j)=0
      Do k=1,4
      B(i,j)=B(i,j)+AMAT(i,k)*GMAT(k,j)
      Enddo
      Enddo
      Enddo
      Do i=1,8
      B(4,2*i-1)=SH(i)/rad
      B(4,2*i)=0
      Enddo
      Call dmat_0()
      Do i=1,4
      Do j=1,16
      DB(i,j)=0
      Do k=1,4
      DB(i,j)=DB(i,j)+D(i,k)*B(k,j)
      Enddo
      Enddo
      Enddo
      If(issol.eq.1)Then
      Do i=1,8
      iin=ndn*(NOC(eind,i)-1)
      ii=ndn*(i-1)
      Do j=1,ndn
      Q(ii+j)=F(iin+j)
      Enddo
      Enddo
      SHT(1)=0.25*(1-xi)*(1-eta)
      SHT(2)=0.25*(1+xi)*(1-eta)
      SHT(3)=0.25*(1+xi)*(1+eta)
      SHT(4)=0.25*(1-xi)*(1+eta)
c---随着单元的变化，温度增量需要变化
      dte=0
      tem=0
      Do i=1,4
      dte=dte+SHT(i)*DT(NODES(i))
      tem=tem+SHT(i)*T00(NODES(i))
      Enddo
c      alpha=PM(6)
      alpha=thermalpha(tem)
      c1=alpha*dte
      Do i=1,4
      c=0
      Do k=1,16
      c=c+DB(i,k)*Q(k)
      Enddo
c---在节点位移已知情况下，代入的是高斯积分点的坐标，求出的应力应该是积分点处的应力此处有4个
      STR(eind,NIK,i)=c-c1*(D(i,1)+D(i,2)+D(i,4))
      Enddo
      Endif
      Return     
      End Subroutine
      
      Subroutine tymath_1(eind,n)
      Use Tsd
      Implicit None
      Integer(kind=4) :: i,j,eind,n
      j2=1.0/6*((STR(eind,n,1)-STR(eind,n,4))**2+(STR(eind,n,1)-
     1STR(eind,n,2))
     1**2+(STR(eind,n,4)-STR(eind,n,2))**2)+STR(eind,n,3)**2
      efs=j2**0.5-strainhard*INDENS(eind,n)**0.5
      if(efs.lt.0)efs=0
      g=efs**strexp/(2*j2**0.5)
      strdev(1)=STR(eind,n,1)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1STR(eind,n,4))
      strdev(2)=STR(eind,n,2)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1STR(eind,n,4))
      strdev(3)=STR(eind,n,3)
      strdev(4)=STR(eind,n,4)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1STR(eind,n,4))
      dj2dsigma(1)=2.0/3*STR(eind,n,1)-1.0/3*STR(eind,n,4)-
     11.0/3*STR(eind,n,2)
      dj2dsigma(2)=2.0/3*STR(eind,n,2)-1.0/3*STR(eind,n,4)-
     11.0/3*STR(eind,n,1)
      dj2dsigma(3)=2.0*STR(eind,n,3)
      dj2dsigma(4)=2.0/3*STR(eind,n,4)-1.0/3*STR(eind,n,1)-
     11.0/3*STR(eind,n,2)
      do i=1,4
      f1(i)=bv0*INDENS(eind,n)*exp(-acten/boltcon/tem)
     1*g*strdev(i)
      Enddo
      f2=kv0*INDENS(eind,n)*exp(-acten/boltcon/tem)*efs**
     1(strexp+matecon)
      End Subroutine

      subroutine tymath_2(eind,n)
c----各种某时刻时某一个单元eind上的某一点NIX(&n)的那些等式及其导数值-----
c----此处应该是为需要计算每个单元的应力函数矩阵等服务的
c----此函数涉及到初始时刻的应力值，需要计算初始时刻的内应力
      use Tsd
      Implicit None
      Integer(kind=4) :: i,j,eind,n
c---函数parainitl()放在这里是由于以后计算降温过程中参数是有变化的
c----	call parainitl()，参数先给的是初值，然后如果考虑变化，
c----则可在此设置一个计算参数变化的函数
      j2=1.0/6*((STR(eind,n,1)-STR(eind,n,4))**2+(STR(eind,n,1)-
     1  STR(eind,n,2))**2+
     1  (STR(eind,n,4)-STR(eind,n,2))**2)+STR(eind,n,3)**2
      efs=j2**0.5-strainhard*INDENS(eind,n)**0.5
      If(efs.lt.0)efs=0
      g=efs**strexp/(2*j2**0.5)
      strdev(1)=STR(eind,n,1)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1  STR(eind,n,4))
      strdev(2)=STR(eind,n,2)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1  STR(eind,n,4))
      strdev(3)=STR(eind,n,3)
      strdev(4)=STR(eind,n,4)-1.0/3*(STR(eind,n,1)+STR(eind,n,2)+
     1  STR(eind,n,4))
      dj2dsigma(1)=2.0/3*STR(eind,n,1)-1.0/3*STR(eind,n,4)-
     1  1.0/3*STR(eind,n,2)
      dj2dsigma(2)=2.0/3*STR(eind,n,2)-1.0/3*STR(eind,n,4)-
     1  1.0/3*STR(eind,n,1)
      dj2dsigma(3)=2.0*STR(eind,n,3)
      dj2dsigma(4)=2.0/3*STR(eind,n,4)-1.0/3*STR(eind,n,1)-
     1  1.0/3*STR(eind,n,2)
      do i=1,4
      f1(i)=bv0*INDENS(eind,n)*exp(-acten/boltcon/tem)
     1  *g*strdev(i)
      enddo
      f2=kv0*INDENS(eind,n)*exp(-acten/boltcon/tem)*efs**
     1  (strexp+matecon)
      dgdj2=efs**(strexp-1)/(4*j2)*(strexp-efs*j2**(-0.5))
      dgdnm=-strainhard*strexp*efs**(strexp-1)/(4*INDENS(eind,n)
     1  **0.5*j2**0.5)
      df1dsigma(1,1)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(1)*dj2dsigma(1)+g*2.0/3)
      df1dsigma(1,2)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(1)*dj2dsigma(2)+g*(-1.0/3))
      df1dsigma(1,3)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(1)*dj2dsigma(3)+g*0.0)
      df1dsigma(1,4)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(1)*dj2dsigma(4)+g*(-1.0/3))
      df1dsigma(2,1)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(2)*dj2dsigma(1)+g*(-1.0/3))
      df1dsigma(2,2)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(2)*dj2dsigma(2)+g*2.0/3)
      df1dsigma(2,3)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(2)*dj2dsigma(3)+g*0.0)
      df1dsigma(2,4)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(2)*dj2dsigma(4)+g*(-1.0/3))
      df1dsigma(3,1)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(3)*dj2dsigma(1)+g*0.0)
      df1dsigma(3,2)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(3)*dj2dsigma(2)+g*0.0)
      df1dsigma(3,3)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(3)*dj2dsigma(3)+g*1.0)
      df1dsigma(3,4)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(3)*dj2dsigma(4)+g*0.0)
      df1dsigma(4,1)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(4)*dj2dsigma(1)+g*(-1.0/3))
      df1dsigma(4,2)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(4)*dj2dsigma(2)+g*(-1.0/3))
      df1dsigma(4,3)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(4)*dj2dsigma(3)+g*0.0)
      df1dsigma(4,4)=bv0*exp(-acten/boltcon/tem)*INDENS(eind,n)*
     1  (dgdj2*strdev(4)*dj2dsigma(4)+g*2.0/3)
      do i=1,4
      df1dnm(i)=bv0*exp(-acten/boltcon/tem)*(g*strdev(i)+
     1  INDENS(eind,n)*dgdnm*strdev(i))
      df1dt(i)=bv0*INDENS(eind,n)*g*exp(-acten/boltcon/tem)*
     1  acten/boltcon/(tem**2)*strdev(i)
      df2dsigma(i)=kv0*INDENS(eind,n)*
     1  exp(-acten/boltcon/tem)*(strexp+matecon)
     1  *efs**(strexp+matecon-1)/(2*j2**0.5)*dj2dsigma(i)
      enddo
      df2dnm=kv0*exp(-acten/boltcon/tem)*efs**
     1  (strexp+matecon-1)*(efs-
     1  INDENS(eind,n)**0.5*(strexp+matecon)*strainhard/2.0)
      df2dt=kv0*INDENS(eind,n)*efs**(strexp+matecon)*
     1  exp(-acten/boltcon/tem)*acten/boltcon/tem**2
      do i=1,4
      do j=1,4
      hn(i,j)=Theta*DeltT/(1-Theta*DeltT*df2dnm)*df1dnm(i)*
     1  df2dsigma(j)+df1dsigma(i,j)
      enddo
c talpha(i):may be change to call thmalpha(tem)
c      talpha(i)=PM(6)
      enddo
      talpha(1)=thermalpha(tem)
      talpha(2)=talpha(1)
      talpha(3)=0
      talpha(4)=talpha(1)
      end subroutine

      Subroutine dbtmat_2(eind)
c---原始的应力函数矩阵
c---此函数是有功能转换的issol
      Use Tsd
      Implicit None
      Integer(kind=4)::eind
      Integer(kind=4):: i,j,k,iin,ii
      Integer(kind=4):: NODES(8)
      Real(kind=8)::dte,ddte
      Real(kind=8)::pi,xi,eta,c,c1,c2,Tm
      Real(kind=8)::SH(8),SHX(8),SHY(8),SHT(4),R(8),Z(8),TMAT(4)
      Real(kind=8)::AMAT(3,4),GMAT(4,16),DDTmat(4),D1(4,4),D1D(4,16)
      pi=3.1415926
      Tm=800
      xi=NIX(NIK,1)
      eta=NIX(NIK,2)
      Do i=1,8
      NODES(i)=NOC(eind,i)
      Enddo
c---R,Z,T
      Do i=1,8
      R(i)=X(NODES(i),1)
      Z(i)=X(NODES(i),2)
      Enddo
      Do i=1,4
      TMAT(i)=T01(NODES(i))
      DDTmat(i)=DDT(NODES(i))
      Enddo
c---shape function-----
      SH(1)=-0.25*(1-xi)*(1-eta)*(1+xi+eta)
      SH(2)=-0.25*(1+xi)*(1-eta)*(1-xi+eta)
      SH(3)=-0.25*(1+xi)*(1+eta)*(1-xi-eta)
      SH(4)=-0.25*(1-xi)*(1+eta)*(1+xi-eta)
      SH(5)=0.5*(1-xi*xi)*(1-eta)
      SH(6)=0.5*(1+xi)*(1-eta*eta)
      SH(7)=0.5*(1-xi*xi)*(1+eta)
      SH(8)=0.5*(1-xi)*(1-eta*eta)
c---shape derivative function
      SHX(1)=0.25*(1-eta)*(2*xi+eta)
      SHY(1)=0.25*(1-xi)*(2*eta+xi)
      SHX(2)=0.25*(1-eta)*(2*xi-eta)
      SHY(2)=0.25*(1+xi)*(2*eta-xi)
      SHX(3)=0.25*(1+eta)*(2*xi+eta)
      SHY(3)=0.25*(1+xi)*(2*eta+xi)
      SHX(4)=0.25*(1+eta)*(2*xi-eta)
      SHY(4)=0.25*(1-xi)*(2*eta-xi)
      SHX(5)=-xi*(1-eta)
      SHY(5)=-0.5*(1-xi*xi)
      SHX(6)=0.5*(1-eta*eta)
      SHY(6)=-eta*(1+xi)
      SHX(7)=-xi*(1+eta)
      SHY(7)=0.5*(1-xi*xi)
      SHX(8)=-0.5*(1-eta*eta)
      SHY(8)=-eta*(1-xi)
c---由于位错中实际温度的存在，4点形函数
      SHT(1)=0.25*(1-xi)*(1-eta)
      SHT(2)=0.25*(1+xi)*(1-eta)
      SHT(3)=0.25*(1+xi)*(1+eta)
      SHT(4)=0.25*(1-xi)*(1+eta)
c---formation of jacobian
      TJ11=0
      TJ12=0
      TJ21=0
      TJ22=0
      Do i=1,8
      TJ11=TJ11+SHX(i)*R(i)
      TJ12=TJ12+SHX(i)*Z(i)
      TJ21=TJ21+SHY(i)*R(i)
      TJ22=TJ22+SHY(i)*Z(i)
      Enddo
c---determinant of jacobian
      DJ=TJ11*TJ22-TJ12*TJ21
c---A,G matrix ,geometry function
      Do i=1,3
      Do j=1,4
      AMAT(i,j)=0
      Enddo
      Enddo
      AMAT(1,1)=TJ22/DJ
      AMAT(1,2)=-TJ12/DJ
      AMAT(2,3)=-TJ21/DJ
      AMAT(2,4)=TJ11/DJ
      AMAT(3,1)=-TJ21/DJ
      AMAT(3,2)=TJ11/DJ
      AMAT(3,3)=TJ22/DJ
      AMAT(3,4)=-TJ12/DJ
      Do i=1,4
      Do j=1,16
      GMAT(i,j)=0
      Enddo
      Enddo
      Do i=1,8
      GMAT(1,2*i-1)=SHX(i)
      GMAT(2,2*i-1)=SHY(i)
      GMAT(3,2*i)=SHX(i)
      GMAT(4,2*i)=SHY(i)
      Enddo
c---任意一点处实际温度tem以及半径rad的表示
      rad=0
      tem=0
      Do i=1,8
      rad=rad+SH(i)*R(i)
      enddo
      Do i=1,4
      tem=tem+SHT(i)*TMAT(i)
      Enddo  
c---B matrix------------------
      Do i=1,3
      Do j=1,16
      B(i,j)=0
      Do k=1,4
      B(i,j)=B(i,j)+AMAT(i,k)*GMAT(k,j)
      Enddo
      Enddo
      Enddo
      Do i=1,8
      B(4,2*i-1)=SH(i)/rad
      B(4,2*i)=0
      Enddo
c*********************~D**************************
      Call tymath_2(eind,NIK)
      Call dmat_0()
      Do i=1,4
      Do j=1,4
      D1(i,j)=0
      Do k=1,4
      D1(i,j)=D1(i,j)+Theta*DeltT*D(i,k)*hn(k,j)
      Enddo
      Enddo
      Enddo

      Do i=1,4
      D1(i,i)=D1(i,i)+1
      Enddo
      
      Do i=1,4
      Do j=1,8
      D1D(i,j)=0
      Enddo
      Enddo

      Do i=1,4
      Do j=1,4
      D1D(i,j)=D1(i,j)
      Enddo
      D1D(i,i+4)=1
      Enddo

      Do i=1,4
      Do j=i,4
      If(D1D(j,i).ne.0)Then
      Do k=i+1,8
      D1D(j,k)=D1D(j,k)/D1D(j,i)
      Enddo
      D1D(j,i)=1
      Endif
      Enddo
      Do j=i+1,4
      If(D1D(j,i).ne.0)Then
      Do k=1,8
      D1D(j,k)=D1D(j,k)-D1D(i,k)
      Enddo
      Endif
      Enddo
      Enddo
c**********
c**********
      Do j=4,2,-1
      Do i=1,j-1
      If(D1D(i,j).ne.0)Then
      Do k=j+1,8
      D1D(i,k)=D1D(i,k)-D1D(i,j)*D1D(j,k)
      Enddo
      D1D(i,j)=0
      Endif
      Enddo
      Enddo
c*************************************************
c      Do i=1,4
c      write(*,'(8es12.4)')(D1D(i,j),j=1,8)
c      Enddo
c*************************************************
      Do i=1,4
      Do j=1,4
      DD(i,j)=0
      Do k=1,4
      DD(i,j)=DD(i,j)+D1D(i,k+4)*D(k,j)
      Enddo
      Enddo
      Enddo
      Do i=1,4
      Do j=1,16
      DB(i,j)=0
      Do k=1,4
      DB(i,j)=DB(i,j)+DD(i,k)*B(k,j)
      Enddo
      Enddo
      Enddo

      If(issol.eq.1)then
      Do i=1,8
      iin=ndn*(NOC(eind,i)-1)
      ii=ndn*(i-1)
      Do j=1,ndn
      Q(ii+j)=F(iin+j)
      Enddo
      Enddo
c-----two types of definition of temperature:dte,ddte
c-----here,c : D*B*Q, total stress part
c ddte,dte:may be changed later
      dte=0
      ddte=0
      Do i=1,4
      ddte=ddte+DDTmat(i)*SHT(i)
      dte=dte+SHT(i)*(Tm-T00(NODES(i)))/Ttol*DeltT
      Enddo
      if(NT2.ne.0)then
      dte=0
      endif
      Do i=1,4
      c=0
      Do k=1,16
      c=c+DB(i,k)*Q(k)
      Enddo
c-----dislocations density model, adapting simulation 1(display)
      c1=0
      c2=(f2*DeltT+Theta*DeltT*dte*df2dt)/(1-Theta*DeltT*df2dnm)
      Do j=1,4
      c1=c1+DD(i,j)*(f1(j)*DeltT+talpha(j)*ddte+Theta*DeltT*df1dnm(j)*c2
     1  +Theta*DeltT*df1dt(j)*dte)
      Enddo
      TSTR(eind,NIK,i)=(c-c1)
      STR(eind,NIK,i)=STR(eind,NIK,i)+TSTR(eind,NIK,i)
      Enddo
c---compute the density of dislocations in integer points
      Do i=1,4
      TINDENS(eind,NIK)=TINDENS(eind,NIK)+df2dsigma(i)*
     1  TSTR(eind,NIK,i)*Theta*DeltT/(1-Theta*DeltT*df2dnm)
      Enddo
      TINDENS(eind,NIK)=TINDENS(eind,NIK)+c2
      INDENS(eind,NIK)=INDENS(eind,NIK)+TINDENS(eind,NIK)
      Endif
      Return     
      End Subroutine
        
      Subroutine temstr(eind)
c---温度载荷计算,是一个时刻下的计算不是计算增量跟下面dtemstr()不同
c---dte的计算需要注意
c---引入content中关于碳化硅晶体膨胀系数---
c---此处是对稳定温度计算的膨胀系数取一个定值，其值在mater.dat中
      Use Tsd
      Implicit None
      Integer(kind=4) :: i,j,k,NODES(4)
      Integer(kind=4)::eind
      Real(kind=8) :: pi
      Real(kind=8) :: xi,eta,SHT(4),alpha,dte
      pi=3.1415926
c      alpha=PM(6)
      dte=0
      xi=NIX(NIK,1)
      eta=NIX(NIK,2)
      SHT(1)=0.25*(1-xi)*(1-eta)
      SHT(2)=0.25*(1+xi)*(1-eta)
      SHT(3)=0.25*(1+xi)*(1+eta)
      SHT(4)=0.25*(1-xi)*(1+eta)
      NODES(1)=NOC(eind,1)
      NODES(2)=NOC(eind,2)
      NODES(3)=NOC(eind,3)
      NODES(4)=NOC(eind,4)
      tem=0
      Do i=1,4
      dte=dte+SHT(i)*DT(NODES(i))
      tem=tem+SHT(i)*T00(NODES(i))
      Enddo
      alpha=thermalpha(tem)
      Do i=1,16
      TL(i)=TL(i)+2*pi*rad*DJ*alpha*dte*(DB(1,i)+DB(2,i)+DB(4,i))
      Enddo
c----here, alpha is const.
      Return
      End Subroutine

      Subroutine temstr_d(eind)
c---temperature load increment compute, new method different from dtemstr
c---alpha referred to <<content>>
c---ddte:veriation of heterogeneous degree
      Use Tsd
      Implicit None
      Integer(kind=4) :: i,j,k,NODES(4)
      Integer(kind=4)::eind
      Real(kind=8) :: pi,Tm
      Real(kind=8) :: xi,eta,SHT(4),dte,ddte,c,c1
      Real(kind=8) :: DDTmat(4)
      pi=3.1415926
      Tm=800
      dte=0
      ddte=0
      xi=NIX(NIK,1)
      eta=NIX(NIK,2)
      SHT(1)=0.25*(1-xi)*(1-eta)
      SHT(2)=0.25*(1+xi)*(1-eta)
      SHT(3)=0.25*(1+xi)*(1+eta)
      SHT(4)=0.25*(1-xi)*(1+eta)
      NODES(1)=NOC(eind,1)
      NODES(2)=NOC(eind,2)
      NODES(3)=NOC(eind,3)
      NODES(4)=NOC(eind,4)
      Do i=1,4
      DDTmat(i)=DDT(NODES(i))
      Enddo
      Do i=1,4
      ddte=ddte+SHT(i)*DDTmat(i)
      dte=dte+SHT(i)*(Tm-T00(NODES(i)))/Ttol*DeltT
      Enddo
      if(NT2.ne.0)then
      dte=0
      endif
      c1=(f2*DeltT+Theta*DeltT*dte*df2dt)/(1-Theta*DeltT*df2dnm)
      Do i=1,16
c      TL(i)=TL(i)+2*pi*rad*DJ*PM(6)*ddte*(DB(1,i)+DB(2,i)+DB(4,i))
      TL(i)=TL(i)+2*pi*rad*DJ*talpha(1)*ddte*(DB(1,i)+DB(2,i)+DB(4,i))
      c=0
      Do j=1,4
      c=c+2*pi*rad*DJ*DB(j,i)*(f1(j)*DeltT+Theta*DeltT*c1*df1dnm(j)+
     1  Theta*DeltT*df1dt(j)*dte)
      Enddo
      TL(i)=TL(i)+c
      Enddo
c      do i=1,16
c          do j=1,4
c              DPL(i)=DPL(i)+2*pi*rad*DJ*DB(j,i)*f1(j)*DeltT
c          enddo
c      enddo
c
c      ctempt=(f2*DeltT+Theta*DeltT*ddte*df2dt)/(1-Theta*DeltT*df2dnm)
c      do i=1,16
c          do j=1,4
c              DPL(i)=PL(i)+2*pi*rad*DJ*DB(j,i)*(f1(j)*DeltT+DeltT*ctempt
c     1*df1dnm(j)+Theta*DeltT*df1dt(j)*ddte)
c          enddo
c      enddo
c----here, alpha is const.
      Return
      End Subroutine
  
      Subroutine elek()
      use Tsd
      Implicit None
      Integer(kind=4) :: i,j,k
      Real(kind=8) :: pi=3.1415926
      Do i=1,16
      Do j=1,16
      Do k=1,4
      SE(i,j)=SE(i,j)+2*pi*rad*B(k,i)*DB(k,j)*DJ
      Enddo
      Enddo
      Enddo
      End subroutine


      subroutine eletoall(eind)
c---firstly compute element stiffness;D:represents actual stiffness matrix.
      use Tsd
      Implicit None
      Integer(kind=4)::eind
      Integer(kind=4) :: i,j,i1,j1,ii,jj,nr,nrt,nc,nct
      Do i1=1,8
      nrt=ndn*(NOC(eind,i1)-1)
      Do ii=1,ndn
      i=ndn*(i1-1)+ii
      nr=nrt+ii
      Do j1=1,8
      nct=ndn*(NOC(eind,j1)-1)
      Do jj=1,ndn
      j=ndn*(j1-1)+jj
      nc=nct+jj-nr+1
      If(nc.gt.0)then
      S(nr,nc)=S(nr,nc)+SE(i,j)
      Endif
      Enddo
      Enddo
      F(nr)=F(nr)+TL(i)
c-----F(:) reprecents the value of righ of equations
c-----In former Subroutine ele2toall(eind),F(nr)=F(nr)+TL(i)+PL(i).
c-----TL(i) consistant of the former TL & PL
c-----Now,i believe the former function is wrong.
      Enddo
      Enddo
      Return
      End subroutine

      subroutine zppostback()
c---输出某时刻的应力及位错密度等信息
c---可选择对某一个时刻的计算数据的提取----
      use Tsd
      Implicit None
      Integer(kind=4) :: i,j,k,num,tmp
      character(len=50)::nfile1,nfile2,nfile3
      character(len=100)::nfile4
      if(NT2.eq.0)then
      write(nfile1,'(f11.3)')Tpass
      write(nfile3,'(f11.3)')DeltT
      write(nfile2,*)'input'
c      nfile4=''//trim(adjustl(nfile1))//'_'//trim(adjustl(nfile3))
c     1//'a.dat'
      nfile4=''//trim(adjustl(nfile1))//'_'
     1//'a.dat'
      nfile1=trim(adjustl(nfile1))//'.dat'
      nfile2=trim(adjustl(nfile2))
      if(mod(NT,1000).eq.0)then
      open(10,file=nfile4)
c      open(11,file=nfile2)
      Write(10,'(5f12.4,i8)')Ttol,Tpass,Tpass2,DeltT,Theta,NT
      write(10,'(a)')'variables="r","z","T0","T1","VMStress","Dens",
     1"str1","str2","str3","str4"'
      write(10,301)imax,jmax
      do j=1,jmax
      do i=1,imax
      num=(j-1)*(3*imax-1)+2*i-1
      write(10,300)X(num,1),X(num,2),T00(num),T01(num),VMSTRESS(num),
     1DENS(num),(STRESS(num,k),k=1,4)
      enddo
      enddo
      write(10,'(a)')'variables=
     1"ins1","ins2","ins3","ins4","Indens'
      Do i=1,ne
      Do j=1,4
      Write(10,303)(STR(i,j,k),k=1,4),INDENS(i,j)
      Enddo
      Enddo
      close(10)
      endif
300   format(2f10.5,8es15.5)
301   format(1x,'zone t= "S_Den",i=',i3,' , j=',i3,' , f= point')
302   format(1x,'zone t= "T"')
303   format(5es15.5)
304   format(1x,'zone t= "InS_Den",i=',i3,' , j=',i3,' , f= point')
400   format(f10.5,2es15.5)
      else
      write(nfile1,'(f11.3)')Tpass2
      write(nfile3,'(f11.3)')DeltT
c      nfile4=''//trim(adjustl(nfile1))//'_'//trim(adjustl(nfile3))
c     1//'r.dat'
      nfile4=''//trim(adjustl(nfile1))//'_'
     1//'r.dat'
      nfile1=trim(adjustl(nfile1))//'.dat'
      if(mod(NT2,100).eq.0)then
      open(10,file=nfile4)
      Write(10,'(5f12.4,i8)')Ttol,Tpass,Tpass2,DeltT,Theta,NT2
      write(10,'(a)')'variables="r","z","T0","T1","VMStress","Dens",
     1"str1","str2","str3","str4"'
      write(10,301)imax,jmax
      do j=1,jmax
      do i=1,imax
      num=(j-1)*(3*imax-1)+2*i-1
      write(10,300)X(num,1),X(num,2),T00(num),T01(num),VMSTRESS(num),
     1DENS(num),(STRESS(num,k),k=1,4)
      enddo
      enddo
      close(10)
      endif
      endif 
      end subroutine
