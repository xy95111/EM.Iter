#######EM algorithm
EM.Iter <-
  function(sdata,n.int = 5,order = 3, r1, max.iter = 1000,cov.rate = 0.1){
    # sdata <- mdata; n.int = 5;order = 3; r1 = 0; max.iter = 500;cov.rate = 0.1
    Xp<-sdata$X
    M<-1###协变量维数
    Kn<-n.int+order
    N<-nrow(sdata)
    ##设置EM算法初值
    b10<-rep(0,M) ##beta1 initial value
    b20<-rep(0,M) ##beta2 initial value
    e0<-1 ##eta initial value
    gamma0<-rep(1,Kn)
    ###均匀节点###
    knots.0 <-seq( from=min(sdata$obs), to=max(sdata$obs), by=(max(sdata$obs)-min(sdata$obs))/(n.int+1) )
    knots<-knots.0[c(-1,-length(knots.0))]##去掉首尾，只剩中间的点##
    bl.Ci<-matrix(0,nrow=N,ncol=Kn) ##400*8的矩阵
    bl.Ci[,]<-iSpline(sdata$obs,degree = 2,knots=knots)
    ###  gamma0的真值用最小二乘得到，i样条关于真值的最小二乘解
    #gamma0 = pmax(as.numeric( ginv(bl.Ci[,]%*%t(bl.Ci[,]))%*%( bl.Ci[,]%*%(0.05*sdata$obs) ) ), 1e-5)
    
    ## 右删失数据C的在风险矩阵YCC
    YCC = outer(sdata$obs, sdata$obs, FUN=">=")*1
    ##alpha初值###
    lambdak<-1/length(unique(mdata$obs))
    alpha0 = rep(0,N)
    for (i in 1:N) {
      if(sdata$Delta[i]==1) alpha0[i]<- lambdak
    }
    
    #alpha0 = 1/as.numeric(t(exp(Xp*b20))%*%YCC)
    rr.mat<-function(x) matrix(x,ncol=N,nrow=Kn,byrow=TRUE)#样条按行展
    ## 回归系数beta1的迭代向量
    bt1em = numeric(max.iter)
    ## 样条系数gamma的迭代矩阵, 每一行是一个迭代结果
    gmem = matrix( 0, max.iter, Kn )
    ## 回归系数beta2的迭代向量
    bt2em = numeric(max.iter)
    ## 小h的迭代矩阵
    ## 小h在每个C处都有值, 只是在u=1的C处有跳, 其他点出为0, 每组小h仍有n个值
    hcem = matrix(0, max.iter, N)
    ## 方差的迭代向量, 每1个元素是1次迭代结果
    vfem = numeric(max.iter)
    ## EM算法迭代初值
    ## beta1的初值
    bt1em[1] = b10  
    ## gamma0初值
    gmem[1,] = gamma0
    ## beta2的初值
    bt2em[1] = b20
    ## 小h的初值
    hcem[1,] = alpha0
    # vf的初值
    vfem[1] = e0
    i<-0
    while( i<max.iter){
      i<-i+1
      ####### 所有参数更新一编###
      #exb<-exp(Xp%*%bt1em[i])
      Lambda1.Ci<-bl.Ci%*%matrix(gmem[i,],ncol=1)
      ## 后面它会做分母, 可能会有0, 修正一下, 后面其实就不需要修正了
      Lambda1.Ci[which(Lambda1.Ci==0)] = 1e-5
      ###中间向量###
      exp.b1x = exp(bt1em[i]*Xp)
      LC.exp = Lambda1.Ci* exp.b1x ###
      ## 右删失部分的中间变量
      exp.b2x = exp(bt2em[i]*Xp)
      HC.C = as.numeric( YCC%*%(sdata$Delta*hcem[i,]) )
      HC.exp = HC.C*exp.b2x
      tha = 1/vfem[i]
      L<-200###蒙特卡洛个数
      
      if(r1>0){
        
        bil = rgamma(L, shape=1/vfem[i], rate=1/vfem[i]) ##bil的值（200个）##
        ##############
        ###中间变量####
        ##############
        r.mat<-function(x) matrix(x,ncol=L,nrow=N,byrow=TRUE)
        exp.v1bv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)*(1+r1*LC.exp%*%bil)^(-1/r1)))
        exp.v1v2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*(1+r1*LC.exp%*%bil)^(-1/r1)))
        Ebi1<-exp.v1bv2/exp.v1v2 #第一种情形bi的条件期望
        exp.v1bbv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)*(1+r1*LC.exp%*%bil)^(-1/r1)))
        Ebi2<-exp.v1bbv2/exp.v1bv2#第二种情形bi的条件期望
        tt1<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)))/L
        #tt2<-(1+1/tha*HC.exp)^(-tha-1)
        exp.v2_v1bv2<-tt1 - exp.v1bv2/L
        tt3<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))))/L
        #tt4<-(1+1/tha*HC.exp)^(-tha)
        exp.v2_v1v2<-tt3 - exp.v1v2/L
        Ebi3<-exp.v2_v1bv2/exp.v2_v1v2#第三种情形bi的条件期望
        tt5<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)))/L
        #tt6<-(1+1/tha)*(1+1/tha*HC.exp)^(-tha-2)
        exp.v2_v1bbv2<-tt5- exp.v1bbv2/L
        Ebi4<-exp.v2_v1bbv2/exp.v2_v1bv2#第四种情形bi的条件期望
        ###bi的条件期望###
        Ebi<-(1-sdata$Delta)*(1-sdata$delta)*Ebi1+(sdata$Delta)*(1-sdata$delta)*Ebi2+(1-sdata$Delta)*(sdata$delta)*Ebi3+(sdata$Delta)*(sdata$delta)*Ebi4
        ##############
        ###Ephibi条件期望####
        ##############
        exp.phiv1bv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)*(1+r1*LC.exp%*%bil)^(-1/r1-1)))
        Ephibi1<-exp.phiv1bv2/exp.v1v2
        exp.phiv1bbv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)*(1+r1*LC.exp%*%bil)^(-1/r1-1)))
        Ephibi2<-exp.phiv1bbv2/exp.v1bv2
        #exp.1_phiv1bv21<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)*(1-(1+r1*LC.exp%*%bil)^(-1/r1-1))))/L
        ##显示解##
        exp.1_phiv1bv2<- tt1-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)*((1+r1*LC.exp%*%bil)^(-1/r1-1))))/L
        Ephibi3<-exp.1_phiv1bv2/exp.v2_v1v2
        # exp.1_phiv1bbv21<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)*(1-(1+r1*LC.exp%*%bil)^(-1/r1-1))))/L
        exp.1_phiv1bbv2<-tt5- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)*((1+r1*LC.exp%*%bil)^(-1/r1-1))))/L
        
        
        
        Ephibi4<-exp.1_phiv1bbv2/exp.v2_v1bv2
        ###phibi的条件期望###
        Ephibi<-(1-sdata$Delta)*(1-sdata$delta)*Ephibi1+(sdata$Delta)*(1-sdata$delta)*Ephibi2+(1-sdata$Delta)*(sdata$delta)*Ephibi3+(sdata$Delta)*(sdata$delta)*Ephibi4
        
        ##############
        ###Ezi条件期望####
        ##############
        Ezi<-0
        Ezi3<-LC.exp*tt1/exp.v2_v1v2
        Ezi4<-LC.exp*tt5/exp.v2_v1bv2
        Ezi<-(1-sdata$Delta)*(sdata$delta)*Ezi3+(sdata$Delta)*(sdata$delta)*Ezi4
        inv.LC = 1/Lambda1.Ci
        inv.LC[which(Lambda1.Ci==0)] = 0
        Ezil = (Ezi*inv.LC)%*%t(gmem[i,])*(bl.Ci[,]) ###N*L矩阵
        
        ##############
        ###Elogbi条件期望####
        ##############
        exp.v1logbv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(log(bil))*(1+r1*LC.exp%*%bil)^(-1/r1)))
        Elogbi1<-exp.v1logbv2/exp.v1v2
        exp.v1blogbv2<- rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil*log(bil))*(1+r1*LC.exp%*%bil)^(-1/r1)))
        Elogbi2<-exp.v1blogbv2/exp.v1bv2
        a1<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(log(bil)))/L )
        #a2<-(digamma(tha)-log(tha+HC.exp))*(1+1/tha*HC.exp)^(-tha)
        Elogbi3<-(a1-exp.v1logbv2/L)/exp.v2_v1v2
        a3<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(log(bil)*bil))/L )
        #a4<-(digamma(tha+1)-log(tha+HC.exp))*(1+1/tha*HC.exp)^(-tha-1)
        Elogbi4<-(a3-exp.v1blogbv2/L)/exp.v2_v1bv2
        Elogbi<-(1-sdata$Delta)*(1-sdata$delta)*Elogbi1+(sdata$Delta)*(1-sdata$delta)*Elogbi2+(1-sdata$Delta)*(sdata$delta)*Elogbi3+(sdata$Delta)*(sdata$delta)*Elogbi4
        
        
      }else{
        
        
        bil = rgamma(L, shape=1/vfem[i], rate=1/vfem[i]) ##bil的值（200个）##
        ##############
        ###中间变量####
        ##############
        r.mat<-function(x) matrix(x,ncol=L,nrow=N,byrow=TRUE)
        exp.v1bv2<- rowSums( as.matrix(exp(-(HC.exp+LC.exp)%*%t(bil)))*r.mat(bil))
        exp.v1v2<- rowSums( as.matrix(exp(-(HC.exp+LC.exp)%*%t(bil))))
        Ebi1<-exp.v1bv2/exp.v1v2 #第一种情形bi的条件期望
        exp.v1bbv2<- rowSums( as.matrix(exp(-(HC.exp+LC.exp)%*%t(bil))*r.mat(bil^2)))
        Ebi2<-exp.v1bbv2/exp.v1bv2#第二种情形bi的条件期望
        tt1<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil)))/L
        #tt2<-(1+1/tha*HC.exp)^(-tha-1)
        exp.v2_v1bv2<-tt1- exp.v1bv2/L
        tt3<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))))/L
        #tt4<-(1+1/tha*HC.exp)^(-tha)
        exp.v2_v1v2<-tt3- exp.v1v2/L
        Ebi3<-exp.v2_v1bv2/exp.v2_v1v2#第三种情形bi的条件期望
        
        tt5<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil^2)))/L
        #tt6<-(1+1/tha)*(1+1/tha*HC.exp)^(-tha-2)
        exp.v2_v1bbv2<-tt5- exp.v1bbv2/L
        Ebi4<-exp.v2_v1bbv2/exp.v2_v1bv2#第四种情形bi的条件期望
        ###bi的条件期望###
        Ebi<-(1-sdata$Delta)*(1-sdata$delta)*Ebi1+(sdata$Delta)*(1-sdata$delta)*Ebi2+(1-sdata$Delta)*(sdata$delta)*Ebi3+(sdata$Delta)*(sdata$delta)*Ebi4
        ##############
        ###Ephibi条件期望####
        ##############
        Ephibi<-Ebi
        
        ##############
        ###Ezi条件期望####
        ##############
        Ezi<-0
        Ezi3<-LC.exp*tt1/exp.v2_v1v2
        Ezi4<-LC.exp*tt5/exp.v2_v1bv2
        Ezi<-(1-sdata$Delta)*(sdata$delta)*Ezi3+(sdata$Delta)*(sdata$delta)*Ezi4
        inv.LC = 1/Lambda1.Ci
        inv.LC[which(Lambda1.Ci==0)] = 0
        Ezil = (Ezi*inv.LC)%*%t(gmem[i,])*(bl.Ci[,]) ###N*L矩阵
        
        ##############
        ###Elogbi条件期望####
        ##############
        exp.v1logbv2<- rowSums( as.matrix(exp(-(HC.exp+LC.exp)%*%t(bil))*r.mat(log(bil))))
        Elogbi1<-exp.v1logbv2/exp.v1v2
        exp.v1blogbv2<- rowSums( as.matrix(exp(-(HC.exp+LC.exp)%*%t(bil))*r.mat(bil*log(bil))))
        Elogbi2<-exp.v1blogbv2/exp.v1bv2
        a1<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(log(bil)))/L )
        #a2<-(digamma(tha)-log(tha+HC.exp))*(1+1/tha*HC.exp)^(-tha)
        Elogbi3<-(a1-exp.v1logbv2/L)/exp.v2_v1v2
        a3<-rowSums( as.matrix(exp(-HC.exp%*%t(bil))*r.mat(bil*log(bil)))/L )
        #a4<-(digamma(tha+1)-log(tha+HC.exp))*(1+1/tha*HC.exp)^(-tha-1)
        Elogbi4<-(a3-exp.v1blogbv2/L)/exp.v2_v1bv2
        Elogbi<-(1-sdata$Delta)*(1-sdata$delta)*Elogbi1+(sdata$Delta)*(1-sdata$delta)*Elogbi2+(1-sdata$Delta)*(sdata$delta)*Elogbi3+(sdata$Delta)*(sdata$delta)*Elogbi4
        
        
        
        
        
      }
      #M-step
      ## 参数估计
      ###### 第1部分, beta10 与 样条系数gamma 的估计 ##############
      Ezil.col = colSums( Ezil)
      ## es.afa 是beta1的函数
      es.afa = function(BT)
      {
        Ezil.col/rowSums( rr.mat(exp(BT*Xp)*Ephibi)*t(bl.Ci) )
      }
      
      ## beta的估计方程
      beta1.score = function(BT)
      {
        ss0 = rowSums( rr.mat(exp(BT*Xp)*Ephibi)*t(bl.Ci) )
        ss1 = rowSums( rr.mat(exp(BT*Xp)*Ephibi*Xp)*t(bl.Ci) )
        sum(t(Ezil)%*%Xp) - sum(Ezil%*%(ss1/ss0))
      }
      
      bt1em[i+1] = nleqslv( bt1em[i], beta1.score )$x
      gmem[i+1,] = es.afa(bt1em[i+1])
      
      ###### 第2部分, beta20 与 alpha0 的估计 ##############
      
      es.hc = function(GM)
      {
        sdata$Delta/ as.numeric(t(exp(GM*Xp)*Ebi)%*%YCC)
      }
      
      
      
      ## 中间变量
      sum.uZ = sum(sdata$Delta*Xp)
      
      ## beta2的估计方程
      beta2.score = function(GM)
      {
        exp.g.Eb = exp(GM*Xp)*Ebi
        ss0 = as.numeric(t(exp.g.Eb)%*%YCC)
        ss1 = as.numeric(t(exp.g.Eb*Xp)%*%YCC)
        
        sum.uZ - sum(sdata$Delta*(ss1/ss0))
      }
      
      bt2em[i+1] = nleqslv( bt2em[i], beta2.score )$x
      
      hcem[i+1,] = es.hc(bt2em[i+1])
      
      ###### 第3部分,方差vf 的估计 ##############
      ## 中间变量
      mean.log.bf = mean(Elogbi) - mean(Ebi)
      ## theta的目标函数
      theta.neg.ll = function(TH)
      {
        ll = TH*log(TH) - lgamma(TH) + TH*mean.log.bf
        -ll
      }
      
      htheta = nlm(theta.neg.ll, 1/vfem[i])$estimate
      htheta<-max(htheta,1e-3)
      vfem[i+1] = 1/htheta
      if(  all(abs(bt1em[i+1] - bt1em[i])<1e-2) 
           & all(abs(bt2em[i+1] - bt2em[i])<1e-2) 
           & all(abs(vfem[i+1] - vfem[i])<1e-2 )
           & all(abs(gmem[i+1,] - gmem[i,])<1e-2)
           & all(abs(hcem[i+1,] - hcem[i,])<1e-2)
           | i+1==max.iter  ) break
      
    }
    #list(b10=bt1em[i+1],b20=bt2em[i+1], e0=vfem[i+1] ,i=i)
    return(list(b10=bt1em[i+1],b20=bt2em[i+1], e0=vfem[i+1],gm0=gmem[i+1,],ap0=hcem[i+1,],iter=i+1,iter=i+1) )   
  }  
