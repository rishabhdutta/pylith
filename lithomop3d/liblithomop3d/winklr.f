c -*- Fortran -*-
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c                             Charles A. Williams
c                       Rensselaer Polytechnic Institute
c                        (C) 2004  All Rights Reserved
c
c  Copyright 2004 Rensselaer Polytechnic Institute.
c  All worldwide rights reserved.  A license to use, copy, modify and
c  distribute this software for non-commercial research purposes only
c  is hereby granted, provided that this copyright notice and
c  accompanying disclaimer is not modified or removed from the software.
c
c  DISCLAIMER:  The software is distributed "AS IS" without any express
c  or implied warranty, including but not limited to, any implied
c  warranties of merchantability or fitness for a particular purpose
c  or any warranty of non-infringement of any current or pending patent
c  rights.  The authors of the software make no representations about
c  the suitability of this software for any particular purpose.  The
c  entire risk as to the quality and performance of the software is with
c  the user.  Should the software prove defective, the user assumes the
c  cost of all necessary servicing, repair or correction.  In
c  particular, neither Rensselaer Polytechnic Institute, nor the authors
c  of the software are liable for any indirect, special, consequential,
c  or incidental damages related to the software, to the maximum extent
c  the law permits.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c
      subroutine winklr(alnz,iwink,wink,histry,nstep,nwink,nhist,nnz,
     & lastep,idout,kto,kw)
c
c...program to implement winkler restoring forces on nodes
c   designated by iwink.  this program adds a constant given
c   by wink to the corresponding diagonal elements of the
c   global stiffness matrix.
c
      include "implicit.inc"
c
c...  subroutine arguments
c
      integer nstep,nwink,nhist,nnz,lastep,idout,kto,kw
      integer iwink(2,nwink)
      double precision alnz(nnz),wink(nwink),histry(nhist,lastep+1)
c
c...  local variables
c
      integer i,mode,k,ihist
c
cdebug      write(6,*) "Hello from winklr_f!"
c
      do i=1,nwink
        mode=iwink(1,i)
        k=iwink(2,i)
        if(mode.eq.1) then
          alnz(k)=alnz(k)+wink(i)
        else
          ihist=-mode
          if(ihist.gt.nhist) then
            if(idout.gt.1) write(kw,1000) nhist,i,k
            write(kto,1000) nhist,i,k
            stop
          end if
          alnz(k)=alnz(k)+wink(i)*histry(ihist,nstep+1)
        end if
      end do
      return
1000  format(//' fatal winkler force error!'//
     & ' attempt to use undefined load history # ',i5,
     & ' for winkler force # ',i5,' eqn. # = ',i7)
      end
c
c version
c $Id: winklr.f,v 1.1 2004/04/14 21:18:30 willic3 Exp $
c
c Generated automatically by Fortran77Mill on Wed May 21 14:15:03 2003
c
c End of file 
