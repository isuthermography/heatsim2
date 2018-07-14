import sys
import copy
import numbers
import sys
import numpy as np # not really necessary


if sys.version_info >= (3,0,0):
    basestring=str
    pass
 


            
class linear_expression(object):
    le_cmdlist=None  # RPN list of (OPERATOR,OptionalValue) tuples
    final=None   # Boolean: Are we finalized? (Finalized means le_cmdlist converted to tuple)
    # The ONLY way to make a non-final expression is to initialize
    # with no arguments

    # Operators
    OP_ADD=1
    OP_MUL=2
    OP_PARAM=3   # Loads OptionalValue onto stack
    OP_VAR=4  # Variable... OptionalValue is (name,coefficient,index0,index1)
    OP_GROUP=5 # Top of value stack is a group and cannot be split up
    OP_DIV=6
    
    def __init__(self,*args):
        
        if len(args)==0:
            self.le_cmdlist=[]
            self.final=False
            return
        
        elif len(args)==1:
            value=args[0]
            if hasattr(value,'le_cmdlist'):
                self.le_cmdlist=copy.deepcopy(value.le_cmdlist)
                self.finalize()
                pass
            elif isinstance(value,numbers.Number):
                self.le_cmdlist=[]
                self.le_cmdlist.append((self.OP_PARAM,value))
                self.finalize()
                pass
            elif isinstance(value,basestring):
                self.le_cmdlist=[]
                self.le_cmdlist.append((self.OP_VAR,(value,1.0,None,None)))
                self.finalize()
                pass
            else:
                raise ValueError("Unknown argument type")
            pass
        else:
            raise ValueError("Too many arguments (%d)" % (len(args)))
        pass
    
    def __add__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.append((self.OP_ADD,None))
        res.finalize()
        return res

    def __radd__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.append((self.OP_ADD,None))
        res.finalize()
        return res

    def __sub__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.append((self.OP_PARAM,-1))
        res.le_cmdlist.append((self.OP_MUL,None))
        res.le_cmdlist.append((self.OP_ADD,None))
        res.finalize()
        return res

    def __neg__(self):
        res=linear_expression()
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.append((self.OP_PARAM,-1))
        res.le_cmdlist.append((self.OP_MUL,None))
        res.finalize()
        return res


    def __rsub__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.append((self.OP_PARAM,-1))
        res.le_cmdlist.append((self.OP_MUL,None))
        res.le_cmdlist.append((self.OP_ADD,None))
        res.finalize()
        return res

    def __getitem__(self,index):
        # index this variable or expression of variables by length 2 tuple index
        res=linear_expression()

        res.le_cmdlist=[]
        for entry in self.le_cmdlist:

            #sys.stderr.write("__getitem__(%s,%s)\n" %(str(self),str(index)))
            #assert(len(self.le_cmdlist)==1 and self.le_cmdlist[0][0]==self.OP_VAR)
            #assert(self.le_cmdlist[0][1][2] is None and self.le_cmdlist[0][1][3] is None)
            if entry[0]==self.OP_VAR and entry[1][2] is None and entry[1][3] is None:
                # must not already be indexed
                
                assert(len(index)==2) # must be two subscripts

                res.le_cmdlist.append((self.OP_VAR,(entry[1][0],entry[1][1],index[0],index[1])))
                pass
            else :
                res.le_cmdlist.append(entry)
                pass
            pass
        
        res.finalize()
        return res

    def __mul__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.append((self.OP_MUL,None))
        res.finalize()
        return res

    def __rmul__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.append((self.OP_MUL,None))
        res.finalize()
        return res

    def __div__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.append((self.OP_DIV,None))
        res.finalize()
        return res
        #raise NotImplementedError("... Use multiplication by reciprocal")

    def __rdiv__(self,param):
        res=linear_expression()
        arg=linear_expression(param)
        res.le_cmdlist.extend(arg.le_cmdlist)
        res.le_cmdlist.extend(self.le_cmdlist)
        res.le_cmdlist.append((self.OP_DIV,None))
        res.finalize()
        return res
        #raise NotImplementedError("... Use multiplication by reciprocal")

    def finalize(self):
        self.final=True
        self.le_cmdlist=tuple(self.le_cmdlist)
        pass
    
    def extractparam(self,pos):
        op=self.le_cmdlist[pos][0]

        # 0 operand operators
        if op==self.OP_PARAM or op==self.OP_VAR:
            return (pos,pos)

        if op==self.OP_GROUP:
            (argstart,argend)=self.extractparam(pos-1)
            return(argstart,pos)
        
        # all other operators are binary (two-operand)

        return (self.extractbinaryparams(pos)[0],pos)
        
        
    def extractbinaryparams(self,pos):
        # extract parameters of a binary operator... ADD or MUL or DIV
        
        op=self.le_cmdlist[pos][0]
        assert(op==self.OP_ADD or op==self.OP_MUL or op==self.OP_DIV)

        argend2=pos-1

        # Collect both arguments
        (argstart2,argend2)=self.extractparam(argend2)

        argend1=argstart2-1
        
        (argstart1,argend1)=self.extractparam(argend1)

        return (argstart1,argend1,argstart2,argend2,pos)
    
    def reduce(self,pos):
        assert(not(self.final))
        # reduce to sum-of-products of variables if possible
        #pos=len(self.le_cmdlist)-1
        shift=0

        op=self.le_cmdlist[pos][0]
        if op==self.OP_PARAM or op==self.OP_VAR:

            # Variable multiplied by 0 -> convert to param of 0
            if op==self.OP_VAR and self.le_cmdlist[pos][1][1]==0.0:
                # print("Variable with zero multiplier")
                self.le_cmdlist[pos]=(self.OP_PARAM,0.0)
                
                return (True,shift,pos,pos)
            

            return (False,0,pos,pos)

        if op==self.OP_GROUP:
            # groups are un-breakable by the reduce function
            # Collect both arguments
            (argstart,argend)=self.extractparam(pos-1)
            return (False,0,argstart,pos)

        # print(str(self))
        # print(argend1)
        # print "Reduce: %d ... (%d,%d)" % (op,self.extractparam(pos)[0],self.extractparam(pos)[1])
        
        # all other operators are binary (two-operand)

        changed=False
        argend2=pos-1

        # Collect both arguments
        # print("reduce %s..." % (str(self)))
        (change,shift2,argstart2,argend2)=self.reduce(argend2)
        # print("... to %s" % (str(self)))
        # print("... to %s" % (self.exprstr()))
        
        # apply shifts and change
        pos+=shift2
        shift+=shift2
        changed|=change
        
        argend1=argstart2-1
        
        # print("reduce %s..." % (str(self)))
        (change,shift1,argstart1,argend1)=self.reduce(argend1)
        # print("... to %s" % (str(self)))
        # print("... to %s" % (self.exprstr()))

        # apply shifts and change
        pos+=shift1
        argstart2+=shift1
        argend2+=shift1
        shift+=shift1
        changed|=change

        # print("argstart1=%d; argend1=%d; argstart2=%d;argend2=%d" % (argstart1,argend1,argstart2,argend2))
        
        # Attempt to reduce

        # ... product of something with zero 
        if op==self.OP_MUL and argend2==argstart2 and self.le_cmdlist[argend2][0]==self.OP_PARAM and self.le_cmdlist[argend2][1]==0.0:
            # print("Product * 0")
            # Remove first argument
            for cnt in range(argend1-argstart1+1):
                del self.le_cmdlist[argstart1]
                shift-=1
                pos-=1
                pass
            # Remove second argument
            del self.le_cmdlist[argstart1]
            shift-=1
            pos-=1
            # reassign self to have value zero
            self.le_cmdlist[argstart1]=(self.OP_PARAM,0.0)
            changed=True
            
            return (changed,shift,argstart1,pos)

        # ... product of zero with something 
        # If this 'if' statement errors out, you are probaby providing 
        # matrix thermal conductivities to the isotropic thermal conductivity model.
        # Try using boundary_conducting_anisotropic instead of boundary_conducting
        if op==self.OP_MUL and argend1==argstart1 and self.le_cmdlist[argend1][0]==self.OP_PARAM and self.le_cmdlist[argend1][1]==0.0:

            # print("Product 0 *")
            # Remove first argument
            del self.le_cmdlist[argstart1]
            shift-=1
            pos-=1
            # Remove second argument
            for cnt in range(argend2-argstart2+1):
                del self.le_cmdlist[argstart1]
                shift-=1
                pos-=1
                pass
            # reassign self to have value zero
            self.le_cmdlist[argstart1]=(self.OP_PARAM,0.0)
            changed=True
            
            return (changed,shift,argstart1,pos)

        # ... product of scalar parameters
        if op==self.OP_MUL and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_PARAM  and self.le_cmdlist[argend2][0]==self.OP_PARAM:
            # print("Product of scalars")
            product=self.le_cmdlist[argend1][1]*self.le_cmdlist[argend2][1]
            del self.le_cmdlist[argstart1]
            del self.le_cmdlist[argstart1]
            shift-=2
            pos-=2
            self.le_cmdlist[argstart1]=(self.OP_PARAM,product)
            changed=True
            
            return (changed,shift,argstart1,pos)

        # ... quotient of scalar parameters
        if op==self.OP_DIV and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_PARAM  and self.le_cmdlist[argend2][0]==self.OP_PARAM:
            # print("Quotient of scalars")
            quotient=float(self.le_cmdlist[argend1][1])/self.le_cmdlist[argend2][1]
            del self.le_cmdlist[argstart1]
            del self.le_cmdlist[argstart1]
            shift-=2
            pos-=2
            self.le_cmdlist[argstart1]=(self.OP_PARAM,quotient)
            changed=True
            
            return (changed,shift,argstart1,pos)

        
        # ... product of scalar*variable
        if op==self.OP_MUL and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_PARAM  and self.le_cmdlist[argend2][0]==self.OP_VAR:
            # print("scalar*variable")
            # Pull coefficient into variable
            coefficient=self.le_cmdlist[argend1][1]*self.le_cmdlist[argend2][1][1]
            varname=self.le_cmdlist[argend2][1][0]
            index0=self.le_cmdlist[argend2][1][2]
            index1=self.le_cmdlist[argend2][1][3]
            
            del self.le_cmdlist[argstart1]
            del self.le_cmdlist[argstart1]
            shift-=2
            pos-=2
            changed=True
            # print(self.le_cmdlist[argstart1])
            self.le_cmdlist[argstart1]=(self.OP_VAR,(varname,coefficient,index0,index1))            
            return (changed,shift,argstart1,pos)

        
        # ... product of variable*scalar
        if op==self.OP_MUL and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_VAR  and self.le_cmdlist[argend2][0]==self.OP_PARAM:
            # print("variable*scalar")
            # Pull coefficient into variable
            coefficient=self.le_cmdlist[argend2][1]*self.le_cmdlist[argend1][1][1]
            varname=self.le_cmdlist[argend1][1][0]
            index0=self.le_cmdlist[argend1][1][2]
            index1=self.le_cmdlist[argend1][1][3]

            del self.le_cmdlist[argstart2]
            del self.le_cmdlist[argstart2]
            shift-=2
            pos-=2
            changed=True
            self.le_cmdlist[argstart1]=(self.OP_VAR,(varname,coefficient,index0,index1))            
            return (changed,shift,argstart1,pos)

        # ... division of variable/scalar
        if op==self.OP_DIV and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_VAR  and self.le_cmdlist[argend2][0]==self.OP_PARAM:
            # print("variable/scalar")
            # Pull coefficient into variable
            coefficient=self.le_cmdlist[argend1][1][1]/self.le_cmdlist[argend2][1]
            varname=self.le_cmdlist[argend1][1][0]
            index0=self.le_cmdlist[argend1][1][2]
            index1=self.le_cmdlist[argend1][1][3]

            del self.le_cmdlist[argstart2]
            del self.le_cmdlist[argstart2]
            shift-=2
            pos-=2
            changed=True
            self.le_cmdlist[argstart1]=(self.OP_VAR,(varname,coefficient,index0,index1))            
            return (changed,shift,argstart1,pos)

        
        # Product involving a sum as the 2nd parameter... distribute it
        if op==self.OP_MUL and self.le_cmdlist[argend2][0]==self.OP_ADD:
            # print("product of sum")
            # print("argstart1=%d; argend1=%d; argstart2=%d;argend2=%d" % (argstart1,argend1,argstart2,argend2))
            firstparam=self.le_cmdlist[argstart1:(argend1+1)]
            # print("firstparam=%d:%d" % (argstart1,argend1))
            (secondparam1start,secondparam1end,secondparam2start,secondparam2end,cmdidx)=self.extractbinaryparams(argend2)
            secondparam1=self.le_cmdlist[secondparam1start:(secondparam1end+1)]
            # print("secondparam1=%d:%d" % (secondparam1start,secondparam1end))
            secondparam2=self.le_cmdlist[secondparam2start:(secondparam2end+1)]
            # print("secondparam2=%d:%d" % (secondparam2start,secondparam2end))
            # remove all of these elements
            for cnt in range(pos+1-argstart1):
                del self.le_cmdlist[argstart1]
                shift-=1
                pos-=1
                pass


            # Distribute sum
            newops=[]
            newops.extend(firstparam)
            newops.extend(secondparam1)
            newops.append((self.OP_MUL,None))
            newops.extend(firstparam)
            newops.extend(secondparam2)
            newops.append((self.OP_MUL,None))
            newops.append((self.OP_ADD,None))

            # insert into cmdlist
            self.le_cmdlist[argstart1:argstart1]=newops
            shift+=len(newops)
            pos+=len(newops)
            changed=True
            return (changed,shift,argstart1,pos)


        # Product involving a sum as the 1st parameter... distribute it
        if op==self.OP_MUL and self.le_cmdlist[argend1][0]==self.OP_ADD:
            # print("sum of product")
            (firstparam1start,firstparam1end,firstparam2start,firstparam2end,cmdidx)=self.extractbinaryparams(argend1)
            firstparam1=self.le_cmdlist[firstparam1start:(firstparam1end+1)]
            firstparam2=self.le_cmdlist[firstparam2start:(firstparam2end+1)]
            secondparam=self.le_cmdlist[argstart2:(argend2+1)]

            # remove all of these elements
            for cnt in range(pos+1-argstart1):
                del self.le_cmdlist[argstart1]
                shift-=1
                pos-=1
                pass


            # Distribute sum
            newops=[]
            newops.extend(firstparam1)
            newops.extend(secondparam)
            newops.append((self.OP_MUL,None))
            newops.extend(firstparam2)
            newops.extend(secondparam)
            newops.append((self.OP_MUL,None))
            newops.append((self.OP_ADD,None))

            # insert into cmdlist
            self.le_cmdlist[argstart1:argstart1]=newops
            shift+=len(newops)
            pos+=len(newops)
            changed=True
            return (changed,shift,argstart1,pos)


        # Quotient involving a sum as the 1st parameter... distribute it
        if op==self.OP_DIV and self.le_cmdlist[argend1][0]==self.OP_ADD:
            # print("sum of product")
            (firstparam1start,firstparam1end,firstparam2start,firstparam2end,cmdidx)=self.extractbinaryparams(argend1)
            firstparam1=self.le_cmdlist[firstparam1start:(firstparam1end+1)]
            firstparam2=self.le_cmdlist[firstparam2start:(firstparam2end+1)]
            secondparam=self.le_cmdlist[argstart2:(argend2+1)]

            # remove all of these elements
            for cnt in range(pos+1-argstart1):
                del self.le_cmdlist[argstart1]
                shift-=1
                pos-=1
                pass


            # Distribute sum
            newops=[]
            newops.extend(firstparam1)
            newops.extend(secondparam)
            newops.append((self.OP_DIV,None))
            newops.extend(firstparam2)
            newops.extend(secondparam)
            newops.append((self.OP_DIV,None))
            newops.append((self.OP_ADD,None))

            # insert into cmdlist
            self.le_cmdlist[argstart1:argstart1]=newops
            shift+=len(newops)
            pos+=len(newops)
            changed=True
            return (changed,shift,argstart1,pos)


        
        # Sum of params... perform addition
        if op==self.OP_ADD and argend1==argstart1 and argend2==argstart2 and self.le_cmdlist[argend1][0]==self.OP_PARAM  and self.le_cmdlist[argend2][0]==self.OP_PARAM:
            # print("Sum of scalars")
            sum=self.le_cmdlist[argend1][1]+self.le_cmdlist[argend2][1]
            del self.le_cmdlist[argstart1]
            del self.le_cmdlist[argstart1]
            shift-=2
            pos-=2
            self.le_cmdlist[argstart1]=(self.OP_PARAM,sum)
            changed=True
            
            return (changed,shift,argstart1,pos)

        
        return (changed,shift,argstart1,pos)


    def fullreduce(self):
        # print("fullreduce()")
        res=linear_expression()
        res.le_cmdlist=list(self.le_cmdlist)
        #oldlist=copy.deepcopy(res.le_cmdlist)
        #print("reducing %s..." % (res.exprstr()))
        (changed,shift,argstart,argend)=res.reduce(len(res.le_cmdlist)-1)
        #print("... to %s" % (str(res)))
        #print("... to %s" % (res.exprstr()))

        while changed:  # res.le_cmdlist != oldlist:
            #print("reducing %s..." % (res.exprstr()))
            #oldlist=copy.deepcopy(res.le_cmdlist)
            (changed,shift,argstart,argend)=res.reduce(len(res.le_cmdlist)-1)
            #if changed: 
            #    #print("... to %s" % (res.exprstr()))
            #    pass
            #else:
            #    print("done.")
            #    pass
            #pass
        res.finalize()
        return res

    def dictform(self):
        res=self.fullreduce()

        mydict={}
        # must first consist of OP_ADD, OP_PARAM, and OP_VAR
        for pos in range(len(res.le_cmdlist)):
            if (res.le_cmdlist[pos][0] != res.OP_ADD and
                res.le_cmdlist[pos][0] != res.OP_PARAM and
                res.le_cmdlist[pos][0] != res.OP_VAR):
                raise ValueError("Expressions with operand %d cannot be reduced to dictionary form" % (res.le_cmdlist[pos][0]))
            if (res.le_cmdlist[pos][0]==res.OP_VAR and (
                    res.le_cmdlist[pos][1][2] is not None or
                    res.le_cmdlist[pos][1][3] is not None)):
                raise ValueError("Expressions with indexed operands such as %s cannot be reduced to dictionary form" % (res.le_cmdlist[pos][1][0]))
            
                    
            pass

        mydict={}
        
        for pos in range(len(res.le_cmdlist)):
            if res.le_cmdlist[pos][0] == res.OP_ADD:
                continue
            elif res.le_cmdlist[pos][0] == res.OP_PARAM:
                # constant parameter
                if "" not in mydict:
                    mydict[""]=0.0
                    pass
                mydict[""]+=res.le_cmdlist[pos][1]
                pass
            
            elif res.le_cmdlist[pos][0] == res.OP_VAR:
                if res.le_cmdlist[pos][1][0] not in mydict:
                    mydict[res.le_cmdlist[pos][1][0]]=0.0
                    pass
                
                mydict[res.le_cmdlist[pos][1][0]]+=res.le_cmdlist[pos][1][1]
                pass
            
            pass
        return mydict
    
    def subexprstr(self,pos):
        if self.le_cmdlist[pos][0]==self.OP_ADD:
            (param2,param2startpos)=self.subexprstr(pos-1)
            (param1,param1startpos)=self.subexprstr(param2startpos-1)
            exprstr="(%s + %s)" % (param1,param2)
            return (exprstr,param1startpos)
        elif self.le_cmdlist[pos][0]==self.OP_MUL:
            (param2,param2startpos)=self.subexprstr(pos-1)
            (param1,param1startpos)=self.subexprstr(param2startpos-1)
            exprstr="(%s * %s)" % (param1,param2)
            return (exprstr,param1startpos)
        elif self.le_cmdlist[pos][0]==self.OP_DIV:
            (param2,param2startpos)=self.subexprstr(pos-1)
            (param1,param1startpos)=self.subexprstr(param2startpos-1)
            exprstr="(%s / %s)" % (param1,param2)
            return (exprstr,param1startpos)
        elif self.le_cmdlist[pos][0]==self.OP_PARAM:
            if isinstance(self.le_cmdlist[pos][1],np.ndarray):
                exprstr="%s" % (str(self.le_cmdlist[pos][1]))
                pass
            else:
                exprstr="%g" % (self.le_cmdlist[pos][1])
                pass
            return (exprstr,pos)
        elif self.le_cmdlist[pos][0]==self.OP_VAR:
            if self.le_cmdlist[pos][1][2] is None and self.le_cmdlist[pos][1][3] is None:
                exprstr="%g*%s" % (self.le_cmdlist[pos][1][1],self.le_cmdlist[pos][1][0])
                pass
            else: 
                exprstr="%g*%s[%d,%d]" % (self.le_cmdlist[pos][1][1],self.le_cmdlist[pos][1][0],self.le_cmdlist[pos][1][2],self.le_cmdlist[pos][1][3]) 
                pass
            return (exprstr,pos)
        elif self.le_cmdlist[pos][0]==self.OP_GROUP:
            (param,paramstartpos)=self.subexprstr(pos-1)
            exprstr="GROUP(%s)" % (param)
            return (exprstr,paramstartpos)
        else:
            raise ValueError("Unknown operator: %d" % (self.le_cmdlist[pos][0]))
        pass
    
    
    def exprstr(self):

        (exprstr,startpos)=self.subexprstr(len(self.le_cmdlist)-1)
        if (startpos!=0):
            raise ValueError("Invalid expression %s. Reduces to %s with %d commands left" % (str(self),exprstr,startpos))
        return exprstr
    

    def __hash__(self):
        # Only finalized expressions can be used as dictionary indexes
        assert(self.final)
        return hash(self.le_cmdlist)

    def __eq__(self,other):
        assert(self.final)
        return self.le_cmdlist == other.le_cmdlist

    def __ne__(self,other):
        assert(self.final)
        return self.le_cmdlist != other.le_cmdlist

        
    def __str__(self):
        res="linear_expression( "
        for pos in range(len(self.le_cmdlist)):
            if self.le_cmdlist[pos][0]==self.OP_PARAM:
                res+="%.8g " % (self.le_cmdlist[pos][1])
                pass
            elif self.le_cmdlist[pos][0]==self.OP_VAR:
                if self.le_cmdlist[pos][1][2] is None and self.le_cmdlist[pos][1][3] is None:
                    res+="%.8g*%s " % (self.le_cmdlist[pos][1][1],self.le_cmdlist[pos][1][0])
                    pass
                else: 
                    res+="%.8g*%s[%d,%d] " % (self.le_cmdlist[pos][1][1],self.le_cmdlist[pos][1][0],self.le_cmdlist[pos][1][2],self.le_cmdlist[pos][1][3]) 
                    pass
                pass
            elif self.le_cmdlist[pos][0]==self.OP_ADD:
                res+="+ "
                pass
            elif self.le_cmdlist[pos][0]==self.OP_MUL:
                res+="* "
                pass
            elif self.le_cmdlist[pos][0]==self.OP_DIV:
                res+="/ "
                pass
            elif self.le_cmdlist[pos][0]==self.OP_GROUP:
                res+="group "
                pass
            else:
                raise ValueError("Unknown operand: %d" % (self.le_cmdlist[pos][0]))
            
            pass
        res+=")"
        return res

        
    pass


def group(linear_exp):
    # Group an expression in a way that is not breakable
    # by reduce() 
    groupval=linear_expression()
    groupval.le_cmdlist=list(linear_exp.le_cmdlist)
    groupval.le_cmdlist.append((groupval.OP_GROUP,None))
    groupval.finalize()
    return groupval

def eliminate_groups(expr):
    groupless=linear_expression()
    groupless.le_cmdlist=[ entry for entry in expr.le_cmdlist if entry[0] != expr.OP_GROUP ]
    groupless.finalize()
    return groupless

def no_groups(expr):
    # return True if expr does not have any groups

    return not(any([ entry[0]==expr.OP_GROUP for entry in expr.le_cmdlist ]))
    

def subst(expr,origvar,newvar_or_value):
    new=linear_expression()
    if isinstance(newvar_or_value,basestring):
        # replace with variable
        new.le_cmdlist=[ (expr.OP_VAR,(newvar_or_value,entry[1][1],entry[1][2],entry[1][3])) if entry[0]==expr.OP_VAR and entry[1][0] == origvar else entry for entry in expr.le_cmdlist ]
        pass
    else:
        # replace with value
        new.le_cmdlist=[]
        for entry in expr.le_cmdlist:
            
            if entry[0]==expr.OP_VAR and entry[1][0]==origvar:
                if entry[1][2] is not None or entry[1][3] is not None:
                    new.le_cmdlist.append((expr.OP_PARAM,newvar_or_value[entry[1][2],entry[1][3]]*entry[1][1]))
                    pass
                else:
                    new.le_cmdlist.append((expr.OP_PARAM,newvar_or_value*entry[1][1]))
                    pass
                pass
            else:
                new.le_cmdlist.append(entry)
                pass
            pass
        pass
    
    new.finalize()
    return new


def crank_subst_in_groups(expr,groupmembers,solnum):
    # Go through expr. Find all groups that contain exclusively
    # variables which are members of the set groupmembers.
    # Substitute (Txxxp + Txxxm)/2 for each (i.e. the Crank-Nicholson
    # replacement. 
    #
    # does not handle nested groups
    # Does not handle variables with indexes

    # print("crank_subst_in_groups...pre: %s\n%s" % (expr.exprstr(),str(expr)))
    groupmembers=set(groupmembers)

    res=linear_expression()
    res.le_cmdlist=list(expr.le_cmdlist)
    
    grouplist=[] # list of (argstart,argend)
    
    pos=len(res.le_cmdlist)-1
    while pos >= 0:
        # print(res.le_cmdlist[pos][0])
        if res.le_cmdlist[pos][0]==res.OP_GROUP:
            # got a group
            # find start of group
            (argstart,argend)=res.extractparam(pos-1)
            grouplist.append((argstart,argend))
            # print("group: %d-%d" % (argstart,argend))
            pos=argstart
            pass
        pos-=1        
        pass

    # Since we started from the end, grouplist is reverse-sorted,
    # with the last entries first
    
    # OK. Take the list of groups and let's see if we should do substitution
    for (argstart,argend) in grouplist:
        #print("Group: %d-%d" % (argstart,argend))
        #print(" ".join([cmd[1][0] for cmd in res.le_cmdlist[argstart:(argend+1)] if cmd[0]==res.OP_VAR]))
        #print("Groupmembers",groupmembers)
        if all([ cmd[1][0] in groupmembers for cmd in res.le_cmdlist[argstart:(argend+1)] if cmd[0]==res.OP_VAR ]):
            #print("Substitute!")
            # all variables in this group are part of groupmembers... do substitution
            # ... since we're doing this from the tail of the expression
            # towards the front (see comment on order of grouplist, above)
            # our insertions don't affect the indices we will be
            # using on the next iteration

            # first remove group command from list
            del res.le_cmdlist[argend+1]

            # Now do substitutions
            for pos in range(argend,argstart-1,-1):
                if res.le_cmdlist[pos][0]==res.OP_VAR:
                    # Found a variable... do substitution!
                    varname=res.le_cmdlist[pos][1][0]
                    coefficient=res.le_cmdlist[pos][1][1]

                    # assume no indexing
                    assert(res.le_cmdlist[pos][1][2] is None and res.le_cmdlist[pos][1][3] is None)

                    # remove old variable
                    del res.le_cmdlist[pos]

                    # insert new  ... average of p and m
                    res.le_cmdlist.insert(pos,(res.OP_VAR,(varname+'p'+str(solnum),coefficient/2.0,None,None)))
                    res.le_cmdlist.insert(pos+1,(res.OP_VAR,(varname+'m',coefficient/2.0,None,None)))
                    res.le_cmdlist.insert(pos+2,(res.OP_ADD,None))
                    pass
                pass
            pass
        pass
    res.finalize()
    # print("crank_subst_in_groups...post: %s\n%s" % (res.exprstr(),str(res)))
    
    return res

if __name__=="__main__":
    # Test linear_expression
    up=linear_expression('up')
    down=linear_expression('down')
    left=linear_expression('left')
    right=linear_expression('right')


    print(str(up*down))

    foo=(up*5-10+down)*30 + (up+5)*(30+12)

    print(str(foo))
    foo=foo.fullreduce()
    print(str(foo))
    
    print(str(foo.dictform()))


    bar=(up*5-group(10+down))*30 + (up+5)*(30+12)
    print(bar.exprstr())
    bar=bar.fullreduce()
    print(str(bar))
    #def heatsim2

    fubar=linear_expression(0)+(linear_expression(6.7)+6.7)*0.5 * -1 * (2*up -4*down + 2*left) + (6.7+6.7)*0.5 * -1 * (2*((up+down)+right) -4*(up-down+(left+right)) + 2*left) + linear_expression('foo!')
    #fubar=(linear_expression(6.0)+6.7)*1 * 1 * (up-1)
    #fubar=(6.7)*0.5 * -1 * (2*up -4*down + 2*left)

    print(str(fubar))
    reduced=fubar.fullreduce()
    print(reduced.exprstr())
    print(str(reduced))
    
    print(str(fubar.dictform()))
    
