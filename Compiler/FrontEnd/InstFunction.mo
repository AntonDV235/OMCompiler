/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

encapsulated package InstFunction
" file:        InstFunction.mo
  package:     InstFunction
  description: Function instantiation

  RCS: $Id: InstFunction.mo 17556 2013-10-05 23:58:57Z adrpo $

  This module is responsible for instantiation of Modelica functions.

"

public import Absyn;
public import ClassInf;
public import Connect;
public import ConnectionGraph;
public import DAE;
public import FCore;
public import InnerOuter;
public import InstTypes;
public import Mod;
public import Prefix;
public import SCode;
public import UnitAbsyn;

protected import Lookup;
protected import MetaUtil;
protected import Inst;
protected import InstBinding;
protected import InstVar;
protected import InstUtil;
protected import UnitAbsynBuilder;
protected import List;
protected import Types;
protected import Flags;
protected import FGraph;
protected import FNode;
protected import Debug;
protected import SCodeDump;
protected import Util;
protected import Config;
protected import DAEUtil;
protected import PrefixUtil;
protected import Error;
protected import Builtin;

protected type Ident = DAE.Ident "an identifier";
protected type InstanceHierarchy = InnerOuter.InstHierarchy "an instance hierarchy";
protected type InstDims = list<list<DAE.Dimension>>;

public function instantiateExternalObject
"instantiate an external object.
 This is done by instantiating the destructor and constructor
 functions and create a DAE element containing these two."
  input FCore.Cache inCache;
  input FCore.Graph inEnv "environment";
  input InnerOuter.InstHierarchy inIH;
  input list<SCode.Element> els "elements";
  input DAE.Mod inMod;
  input Boolean impl;
  input SCode.Comment comment;
  input Absyn.Info info;
  output FCore.Cache outCache;
  output FCore.Graph outEnv;
  output InnerOuter.InstHierarchy outIH;
  output DAE.DAElist dae "resulting dae";
  output ClassInf.State ciState;
algorithm
  (outCache,outEnv,outIH,dae,ciState) := matchcontinue(inCache,inEnv,inIH,els,inMod,impl,comment,info)
    local
      SCode.Element destr,constr;
      FCore.Graph env1;
      FCore.Cache cache;
      Ident className;
      Absyn.Path classNameFQ;
      DAE.Type functp;
      FCore.Graph fs,fs1,env;
      FCore.Ref r;
      InstanceHierarchy ih;
      DAE.ElementSource source "the origin of the element";
      // Explicit instantiation, generate constructor and destructor and the function type.
    case  (cache,env,ih,_,_,false,_,_)
      equation
        className = FNode.refName(FGraph.lastScopeRef(env)); // The external object classname is in top frame of environment.
        checkExternalObjectMod(inMod, className);
        destr = SCode.getExternalObjectDestructor(els);
        constr = SCode.getExternalObjectConstructor(els);
        env = FGraph.mkClassNode(env, destr, Prefix.NOPRE(), inMod);
        env = FGraph.mkClassNode(env, constr, Prefix.NOPRE(), inMod);
        (cache,ih) = instantiateExternalObjectDestructor(cache,env,ih,destr);
        (cache,ih,functp) = instantiateExternalObjectConstructor(cache,env,ih,constr);
        SOME(classNameFQ)=  FGraph.getScopePath(env); // Fully qualified classname
        // Extend the frame with the type, one frame up at the same place as the class.
        (env, r) = FGraph.stripLastScopeRef(env);
        env = FGraph.mkTypeNode(env,className,functp);
        env = FGraph.pushScopeRef(env, r);

        // set the  of this element
       source = DAEUtil.addElementSourcePartOfOpt(DAE.emptyElementSource, FGraph.getScopePath(env));
       source = DAEUtil.addCommentToSource(source, SOME(comment));
       source = DAEUtil.addElementSourceFileInfo(source, info);
      then
        (cache,env,ih,DAE.DAE({DAE.EXTOBJECTCLASS(classNameFQ,source)}),ClassInf.EXTERNAL_OBJ(classNameFQ));

    // Implicit, do not instantiate constructor and destructor.
    case (cache,_,ih,_,_,true,_,_)
      equation
        SOME(classNameFQ)= FGraph.getScopePath(inEnv); // Fully qualified classname
      then
        (cache,inEnv,ih,DAE.emptyDae,ClassInf.EXTERNAL_OBJ(classNameFQ));

    // failed
    else
      equation
        Debug.fprintln(Flags.FAILTRACE, "- InstFunction.instantiateExternalObject failed.");
      then fail();
  end matchcontinue;
end instantiateExternalObject;

protected function checkExternalObjectMod
  "Checks that an external object instance does not have any modifiers. This is
   done because an external object may only have two elements, a constructor and
   a destructor, and there's no point in modifying these."
  input DAE.Mod inMod;
  input String inClassName;
algorithm
  _ := match(inMod, inClassName)
    local
      DAE.Ident id;
      DAE.Mod mod;
      Absyn.Info info;

    case (DAE.NOMOD(), _) then ();
    case (DAE.MOD(subModLst = {}), _) then ();

    // The modifier contains a list of submods. Print an error for the first one
    // to make it look like a normal modifier error.
    case (DAE.MOD(subModLst = DAE.NAMEMOD(ident = id, mod = mod) :: _), _)
      equation
        info = Mod.getModInfo(mod);
        Error.addSourceMessage(Error.MISSING_MODIFIED_ELEMENT,
          {id, inClassName}, info);
      then
        fail();

  end match;
end checkExternalObjectMod;

protected function instantiateExternalObjectDestructor
"instantiates the destructor function of an external object"
  input FCore.Cache inCache;
  input FCore.Graph env;
  input InnerOuter.InstHierarchy inIH;
  input SCode.Element cl;
  output FCore.Cache outCache;
  output InnerOuter.InstHierarchy outIH;
algorithm
  (outCache,outIH) := matchcontinue (inCache,env,inIH,cl)
    local
      FCore.Cache cache;
      FCore.Graph env1;
      InstanceHierarchy ih;

    case (cache,_,ih,_)
      equation
        (cache,_,ih) = implicitFunctionInstantiation(cache,env,ih,DAE.NOMOD(),Prefix.NOPRE(),cl,{});
      then
        (cache,ih);
    else
      equation
        Debug.fprintln(Flags.FAILTRACE, "- InstFunction.instantiateExternalObjectDestructor failed.");
      then fail();
   end matchcontinue;
end instantiateExternalObjectDestructor;

protected function instantiateExternalObjectConstructor
"instantiates the constructor function of an external object"
  input FCore.Cache inCache;
  input FCore.Graph env;
  input InnerOuter.InstHierarchy inIH;
  input SCode.Element cl;
  output FCore.Cache outCache;
  output InnerOuter.InstHierarchy outIH;
  output DAE.Type outType;
algorithm
  (outCache,outIH,outType) := matchcontinue (inCache,env,inIH,cl)
    local
      FCore.Cache cache;
      FCore.Graph env1;
      DAE.Type ty;
      InstanceHierarchy ih;

    case (cache,_,ih,_)
      equation
        (cache,env1,ih) = implicitFunctionInstantiation(cache,env,ih, DAE.NOMOD(), Prefix.NOPRE(), cl, {});
        (cache,ty,_) = Lookup.lookupType(cache,env1,Absyn.IDENT("constructor"),NONE());
      then
        (cache,ih,ty);
    else
      equation
        Debug.fprintln(Flags.FAILTRACE, "- InstFunction.instantiateExternalObjectConstructor failed.");
      then fail();
  end matchcontinue;
end instantiateExternalObjectConstructor;

public function implicitFunctionInstantiation
"This function instantiates a function, which is performed *implicitly*
  since the variables of a function should not be instantiated as for an
  ordinary class."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input DAE.Mod inMod;
  input Prefix.Prefix inPrefix;
  input SCode.Element inClass;
  input list<list<DAE.Dimension>> inInstDims;
  output FCore.Cache outCache;
  output FCore.Graph outEnv;
  output InnerOuter.InstHierarchy outIH;
algorithm
  (outCache,outEnv,outIH):= match (inCache,inEnv,inIH,inMod,inPrefix,inClass,inInstDims)
    local
      DAE.Type ty1;
      FCore.Graph env,cenv;
      Absyn.Path fpath;
      DAE.Mod mod;
      Prefix.Prefix pre;
      SCode.Element c;
      String n;
      InstDims inst_dims;
      FCore.Cache cache;
      InstanceHierarchy ih;
      DAE.ElementSource source "the origin of the element";
      list<DAE.Function> funs;
      DAE.Function fun;
      SCode.Restriction r;
      SCode.Partial pPrefix;

    case (cache,env,ih,mod,pre,(c as SCode.CLASS(name = n,restriction = SCode.R_RECORD(_), partialPrefix = pPrefix)),inst_dims)
      equation
        (cache,c,cenv) = Lookup.lookupRecordConstructorClass(cache,env,Absyn.IDENT(n));
        (cache,env,ih,{DAE.FUNCTION(fpath,_,ty1,_,_,_,_,source,_)}) = implicitFunctionInstantiation2(cache,cenv,ih,mod,pre,c,inst_dims,true);
        // fpath = Absyn.makeFullyQualified(fpath);
        fun = DAE.RECORD_CONSTRUCTOR(fpath,ty1,source,DAE.VARIABLE());
        cache = InstUtil.addFunctionsToDAE(cache, {fun}, pPrefix);
      then (cache,env,ih);

    case (cache,env,ih,mod,pre,(c as SCode.CLASS(restriction = r,partialPrefix = pPrefix)),inst_dims)
      equation
        failure(SCode.R_RECORD(_) = r);
        true = MetaUtil.strictRMLCheck(Flags.isSet(Flags.STRICT_RML),c);
        (cache,env,ih,funs) = implicitFunctionInstantiation2(cache,env,ih,mod,pre,c,inst_dims,false);
        cache = InstUtil.addFunctionsToDAE(cache, funs, pPrefix);
      then (cache,env,ih);

    // handle failure
    case (_,env,_,_,_,SCode.CLASS(name=n),_)
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("- Inst.implicitFunctionInstantiation failed " +& n);
        Debug.traceln("  Scope: " +& FGraph.printGraphPathStr(env));
      then fail();
  end match;
end implicitFunctionInstantiation;

protected function implicitFunctionInstantiation2
"This function instantiates a function, which is performed *implicitly*
  since the variables of a function should not be instantiated as for an
  ordinary class."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input DAE.Mod inMod;
  input Prefix.Prefix inPrefix;
  input SCode.Element inClass;
  input list<list<DAE.Dimension>> inInstDims;
  input Boolean instFunctionTypeOnly "if true, do no additional checking of the function";
  output FCore.Cache outCache;
  output FCore.Graph outEnv;
  output InnerOuter.InstHierarchy outIH;
  output list<DAE.Function> funcs;
algorithm
  (outCache,outEnv,outIH,funcs):= matchcontinue (inCache,inEnv,inIH,inMod,inPrefix,inClass,inInstDims,instFunctionTypeOnly)
    local
      DAE.Type ty,ty1;
      ClassInf.State st;
      FCore.Graph env_1,env,tempenv,cenv;
      Absyn.Path fpath;
      DAE.Mod mod;
      Prefix.Prefix pre;
      SCode.Element c;
      String n;
      InstDims inst_dims;
      SCode.Visibility vis;
      SCode.Partial partialPrefix;
      SCode.Encapsulated encapsulatedPrefix;
      DAE.ExternalDecl extdecl;
      SCode.Restriction restr;
      SCode.ClassDef parts;
      list<SCode.Element> els;
      list<Absyn.Path> funcnames;
      FCore.Cache cache;
      InstanceHierarchy ih;
      DAE.ElementSource source "the origin of the element";
      list<DAE.Element> daeElts;
      list<DAE.Function> resfns;
      list<DAE.FunctionDefinition> derFuncs;
      Absyn.Info info;
      DAE.InlineType inlineType;
      SCode.ClassDef cd;
      Boolean partialPrefixBool, isImpure;
      SCode.Comment cmt;
      SCode.FunctionRestriction funcRest;
      InstTypes.CallingScope cs;
      SCode.Visibility visibility;

    // normal functions
    case (cache,env,ih,mod,pre,SCode.CLASS(classDef=cd, prefixes=SCode.PREFIXES(visibility=visibility), partialPrefix = partialPrefix, name = n,restriction = SCode.R_FUNCTION(funcRest),info = info,cmt=cmt),inst_dims,_)
      equation
        false = SCode.isExternalFunctionRestriction(funcRest);
        isImpure = SCode.isImpureFunctionRestriction(funcRest);

        // if we're not MetaModelica set it to non-partial
        c = Util.if_(Config.acceptMetaModelicaGrammar(),
                     inClass,
                     SCode.setClassPartialPrefix(SCode.NOT_PARTIAL(), inClass));
        cs = Util.if_(instFunctionTypeOnly, InstTypes.TYPE_CALL(), InstTypes.INNER_CALL());
        //print("1 Prefix: " +& PrefixUtil.printPrefixStr(pre) +& " path: " +& n +& "\n");
        (cache,cenv,ih,_,DAE.DAE(daeElts),_,ty,_,_,_) =
          Inst.instClass(cache, env, ih, UnitAbsynBuilder.emptyInstStore(), mod, pre,
            c, inst_dims, true, cs, ConnectionGraph.EMPTY, Connect.emptySet);
        List.map2_0(daeElts,InstUtil.checkFunctionElement,false,info);
        // do not add the stripped class to the env, is already there, not stripped!
        env_1 = env; // Env.extendFrameC(env,c);
        (cache,fpath) = Inst.makeFullyQualified(cache, env_1, Absyn.IDENT(n));
        //print("2 Prefix: " +& PrefixUtil.printPrefixStr(pre) +& " path: " +& Absyn.pathString(fpath) +& "\n");
        cmt = InstUtil.extractClassDefComment(cache, env, cd, cmt);
        derFuncs = InstUtil.getDeriveAnnotation(cd, cmt,fpath,cache,cenv,ih,pre,info);

        (cache) = instantiateDerivativeFuncs(cache,env,ih,derFuncs,fpath,info);

        ty1 = InstUtil.setFullyQualifiedTypename(ty,fpath);
        checkExtObjOutput(ty1,info);
        env_1 = FGraph.mkTypeNode(env_1, n, ty1);

        // set the source of this element
        source = DAEUtil.createElementSource(info, FGraph.getScopePath(env), PrefixUtil.prefixToCrefOpt(pre), NONE(), NONE());
        inlineType = InstUtil.isInlineFunc(c);
        partialPrefixBool = SCode.partialBool(partialPrefix);

        daeElts = InstUtil.optimizeFunctionCheckForLocals(fpath,daeElts,NONE(),{},{},{});
        InstUtil.checkFunctionDefUse(daeElts,info);
        /* Not working 100% yet... Also, a lot of code has unused inputs :( */
        Debug.bcall3(
          false and Config.acceptMetaModelicaGrammar() and not instFunctionTypeOnly,
          InstUtil.checkFunctionInputUsed,daeElts,NONE(),Absyn.pathString(fpath));
      then
        (cache,env_1,ih,{DAE.FUNCTION(fpath,DAE.FUNCTION_DEF(daeElts)::derFuncs,ty1,visibility,partialPrefixBool,isImpure,inlineType,source,SOME(cmt))});

    // External functions should also have their type in env, but no dae.
    case (cache,env,ih,mod,pre,(c as SCode.CLASS(partialPrefix=partialPrefix, prefixes=SCode.PREFIXES(visibility=visibility), name = n,restriction = (restr as SCode.R_FUNCTION(SCode.FR_EXTERNAL_FUNCTION(isImpure))),
        classDef = cd as (parts as SCode.PARTS(elementLst = _)), cmt=cmt, info=info, encapsulatedPrefix = encapsulatedPrefix)),inst_dims,_)
      equation
        (cache,cenv,ih,_,DAE.DAE(daeElts),_,ty,_,_,_) =
          Inst.instClass(cache,env,ih, UnitAbsynBuilder.emptyInstStore(),mod, pre,
            c, inst_dims, true, InstTypes.INNER_CALL(), ConnectionGraph.EMPTY, Connect.emptySet);
        List.map2_0(daeElts,InstUtil.checkFunctionElement,true,info);
        //env_11 = FGraph.mkClassNode(cenv,pre,mod,c);
        // Only created to be able to get FQ path.
        (cache,fpath) = Inst.makeFullyQualified(cache,env,Absyn.IDENT(n));

        cmt = InstUtil.extractClassDefComment(cache, env, cd, cmt);
        derFuncs = InstUtil.getDeriveAnnotation(cd,cmt,fpath,cache,env,ih,pre,info);

        (cache) = instantiateDerivativeFuncs(cache,env,ih,derFuncs,fpath,info);

        ty1 = InstUtil.setFullyQualifiedTypename(ty,fpath);
        checkExtObjOutput(ty1,info);
        ((ty1,_)) = Types.traverseType((ty1,-1),Types.makeExpDimensionsUnknown);
        env_1 = FGraph.mkTypeNode(cenv, n, ty1);
        vis = SCode.PUBLIC();
        (cache,tempenv,ih,_,_,_,_,_,_,_,_,_) =
          Inst.instClassdef(cache, env_1, ih, UnitAbsyn.noStore, mod, pre,
            ClassInf.FUNCTION(fpath,isImpure), n,parts, restr, vis, partialPrefix,
            encapsulatedPrefix, inst_dims, true, InstTypes.INNER_CALL(),
            ConnectionGraph.EMPTY, Connect.emptySet, NONE(), cmt, info) "how to get this? impl" ;
        (cache,ih,extdecl) = instExtDecl(cache, tempenv,ih, n, parts, true, pre,info) "impl" ;

        // set the source of this element
        source = DAEUtil.createElementSource(info, FGraph.getScopePath(env), PrefixUtil.prefixToCrefOpt(pre), NONE(), NONE());
        partialPrefixBool = SCode.partialBool(partialPrefix);
        InstUtil.checkExternalFunction(daeElts,extdecl,Absyn.pathString(fpath));
      then
        (cache,env_1,ih,{DAE.FUNCTION(fpath,DAE.FUNCTION_EXT(daeElts,extdecl)::derFuncs,ty1,visibility,partialPrefixBool,isImpure,DAE.NO_INLINE(),source,SOME(cmt))});

    // Instantiate overloaded functions
    case (cache,env,ih,_,pre,(SCode.CLASS(name = n, prefixes=SCode.PREFIXES(visibility=visibility), restriction = (SCode.R_FUNCTION(SCode.FR_NORMAL_FUNCTION(isImpure))),
          classDef = SCode.OVERLOAD(pathLst = funcnames),cmt=cmt)),_,_)
      equation
        (cache,env,ih,resfns) = instOverloadedFunctions(cache,env,ih,pre,funcnames) "Overloaded functions" ;
        (cache,fpath) = Inst.makeFullyQualified(cache,env,Absyn.IDENT(n));
        resfns = DAE.FUNCTION(fpath,{DAE.FUNCTION_DEF({})},DAE.T_UNKNOWN_DEFAULT,visibility,true,isImpure,DAE.NO_INLINE(),DAE.emptyElementSource,SOME(cmt))::resfns;
      then
        (cache,env,ih,resfns);

    // handle failure
    case (_,env,_,_,_,SCode.CLASS(name=n),_,_)
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("- Inst.implicitFunctionInstantiation2 failed " +& n);
        Debug.traceln("  Scope: " +& FGraph.printGraphPathStr(env));
      then fail();
  end matchcontinue;
end implicitFunctionInstantiation2;

protected function instantiateDerivativeFuncs "instantiates all functions found in derivative annotations so they are also added to the
dae and can be generated code for in case they are required"
  input FCore.Cache cache;
  input FCore.Graph env;
  input InnerOuter.InstHierarchy ih;
  input list<DAE.FunctionDefinition> funcs;
  input Absyn.Path path "the function name itself, must be added to derivative functions mapping to be able to search upwards";
  input Absyn.Info info;
  output FCore.Cache outCache;
algorithm
 // print("instantiate deriative functions for "+&Absyn.pathString(path)+&"\n");
 (outCache) := instantiateDerivativeFuncs2(cache,env,ih,DAEUtil.getDerivativePaths(funcs),path,info);
 // print("instantiated derivative functions for "+&Absyn.pathString(path)+&"\n");
end instantiateDerivativeFuncs;

protected function instantiateDerivativeFuncs2 "help function"
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input list<Absyn.Path> inPaths;
  input Absyn.Path path "the function name itself, must be added to derivative functions mapping to be able to search upwards";
  input Absyn.Info info;
  output FCore.Cache outCache;
algorithm
  (outCache) := matchcontinue(inCache,inEnv,inIH,inPaths,path,info)
    local
      list<DAE.Function> funcs;
      Absyn.Path p;
      FCore.Cache cache;
      FCore.Graph cenv,env;
      InstanceHierarchy ih;
      SCode.Element cdef;
      list<Absyn.Path> paths;
      String fun,scope;

    case(cache,_,_,{},_,_) then (cache);

    // Skipped recursive calls (by looking in cache)
    case(cache,env,ih,p::paths,_,_)
      equation
        (cache,_,cenv) = Lookup.lookupClass(cache,env,p,true);
        (cache,p) = Inst.makeFullyQualified(cache,cenv,p);
        FCore.checkCachedInstFuncGuard(cache,p);
      then instantiateDerivativeFuncs2(cache,env,ih,paths,path,info);

    case(cache,env,ih,p::paths,_,_)
      equation
        (cache,cdef,cenv) = Lookup.lookupClass(cache,env,p,true);
        (cache,p) = Inst.makeFullyQualified(cache,cenv,p);
        // add to cache before instantiating, to break recursion for recursive definitions.
        cache = FCore.addCachedInstFuncGuard(cache,p);
        (cache,_,ih,funcs) =
        implicitFunctionInstantiation2(cache,cenv,ih,DAE.NOMOD(),Prefix.NOPRE(),cdef,{},false);

        funcs = InstUtil.addNameToDerivativeMapping(funcs,path);
        cache = FCore.addDaeFunction(cache, funcs);
      then instantiateDerivativeFuncs2(cache,env,ih,paths,path,info);

    else
      equation
        p :: _ = inPaths;
        fun = Absyn.pathString(p);
        scope = FGraph.printGraphPathStr(inEnv);
        Error.addSourceMessage(Error.LOOKUP_FUNCTION_ERROR,{fun,scope},info);
      then fail();

  end matchcontinue;
end instantiateDerivativeFuncs2;

public function implicitFunctionTypeInstantiation
"author: PA
  When looking up a function type it is sufficient to only instantiate the input and output arguments of the function.
  The implicitFunctionInstantiation function will instantiate the function body, resulting in a DAE for the body.
  This function does not do that. Therefore this function is the only solution available for recursive functions,
  where the function body contain a call to the function itself.

  Extended 2007-06-29, BZ
  Now this function also handles Derived function."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input SCode.Element inClass;
  output FCore.Cache outCache;
  output FCore.Graph outEnv;
  output InnerOuter.InstHierarchy outIH;
algorithm
  (outCache,outEnv,outIH) := matchcontinue (inCache,inEnv,inIH,inClass)
    local
      SCode.Element stripped_class;
      FCore.Graph env_1,env;
      String id,cn2;
      SCode.Partial p;
      SCode.Encapsulated e;
      SCode.Restriction r;
      Option<SCode.ExternalDecl> extDecl;
      list<SCode.Element> elts, stripped_elts;
      FCore.Cache cache;
      InstanceHierarchy ih;
      list<SCode.Annotation> annotationLst;
      Absyn.Info info;
      DAE.DAElist dae;
      list<DAE.Function> funs;
      Absyn.Path cn,fpath;
      Option<list<Absyn.Subscript>> ad;
      SCode.Mod mod1;
      DAE.Mod mod2;
      FCore.Graph cenv;
      SCode.Element c;
      DAE.Type ty1,ty;
      SCode.Prefixes prefixes;
      SCode.Comment cmt;
      list<Absyn.Path> paths;

    // For external functions, include everything essential
    case (cache,env,ih,SCode.CLASS(name = _,
                                   encapsulatedPrefix = _,partialPrefix = _,restriction = SCode.R_FUNCTION(SCode.FR_EXTERNAL_FUNCTION(_)),
                                   classDef = SCode.PARTS(elementLst = _,externalDecl=_)))
      equation
        // stripped_class = SCode.CLASS(id,prefixes,e,p,r,SCode.PARTS(elts,{},{},{},{},{},{},extDecl),cmt,info);
        (cache,env_1,ih,funs) = implicitFunctionInstantiation2(cache, env, ih, DAE.NOMOD(), Prefix.NOPRE(), inClass, {}, true);
        // Only external functions are valid without an algorithm section...
        cache = FCore.addDaeExtFunction(cache, funs);
      then
        (cache,env_1,ih);

    // The function type can be determined without the body. Annotations need to be preserved though.
    case (cache,env,ih,SCode.CLASS(name = id,prefixes = prefixes,
                                   encapsulatedPrefix = e,partialPrefix = p,restriction = _,
                                   classDef = SCode.PARTS(elementLst = elts,externalDecl=_),cmt=cmt, info = info))
      equation
        elts = List.select(elts,isElementImportantForFunction);
        stripped_class = SCode.CLASS(id,prefixes,e,p,SCode.R_FUNCTION(SCode.FR_NORMAL_FUNCTION(false)),SCode.PARTS(elts,{},{},{},{},{},{},NONE()),cmt,info);
        (cache,env_1,ih,_) = implicitFunctionInstantiation2(cache, env, ih, DAE.NOMOD(), Prefix.NOPRE(), stripped_class, {}, true);
        // Only external functions are valid without an algorithm section...
        // cache = FCore.addDaeExtFunction(cache, funs);
      then
        (cache,env_1,ih);

    // Short class definitions.
    case (cache,env,ih,SCode.CLASS(name = id,partialPrefix = _,encapsulatedPrefix = _,restriction = _,
                                   classDef = SCode.DERIVED(typeSpec = Absyn.TPATH(path = cn,arrayDim = _),
                                                            modifications = mod1),info = info))
      equation
        (cache,(c as SCode.CLASS(name = _, restriction = _)),cenv) = Lookup.lookupClass(cache, env, cn, false); // Makes MultiBody gravityacceleration hacks shit itself
        (cache,mod2) = Mod.elabMod(cache, env, ih, Prefix.NOPRE(), mod1, false, Mod.DERIVED(cn), info);

        (cache,_,ih,_,_,_,ty,_,_,_) =
          Inst.instClass(cache,cenv,ih,UnitAbsynBuilder.emptyInstStore(), mod2,
            Prefix.NOPRE(), c, {}, true, InstTypes.INNER_CALL(), ConnectionGraph.EMPTY, Connect.emptySet);

        env_1 = env; // why would you want to do this: FGraph.mkClassNode(env,c); ?????
        (cache,fpath) = Inst.makeFullyQualified(cache,env_1, Absyn.IDENT(id));
        ty1 = InstUtil.setFullyQualifiedTypename(ty,fpath);
        env_1 = FGraph.mkTypeNode(env_1, id, ty1);
        // (cache,env_1,ih,_) = implicitFunctionInstantiation2(cache, env, ih, DAE.NOMOD(), Prefix.NOPRE(), inClass, {}, true);
      then
        (cache,env_1,ih);

    case (cache,env,ih,SCode.CLASS(name = _,partialPrefix = _,encapsulatedPrefix = _,restriction = _,
                                   classDef = SCode.OVERLOAD(pathLst=_)))
      equation
         (cache,env,ih,_) = implicitFunctionInstantiation2(cache, env, ih, DAE.NOMOD(), Prefix.NOPRE(), inClass, {}, true);
      then
        (cache,env,ih);

    case (_,_,_,SCode.CLASS(name=id))
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("- Inst.implicitFunctionTypeInstantiation failed " +& id +& "\nenv: " +& FGraph.getGraphNameStr(inEnv) +& "\nelelement: " +& SCodeDump.unparseElementStr(inClass,SCodeDump.defaultOptions));
      then fail();
  end matchcontinue;
end implicitFunctionTypeInstantiation;

protected function instOverloadedFunctions
"This function instantiates the functions in the overload list of a
  overloading function definition and register the function types using
  the overloaded name. It also creates dae elements for the functions."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input Prefix.Prefix pre;
  input list<Absyn.Path> inAbsynPathLst;
  output FCore.Cache outCache;
  output FCore.Graph outEnv;
  output InnerOuter.InstHierarchy outIH;
  output list<DAE.Function> outFns;
algorithm
  (outCache,outEnv,outIH,outFns) := matchcontinue (inCache,inEnv,inIH,pre,inAbsynPathLst)
    local
      FCore.Graph env,cenv;
      SCode.Element c;
      String id;
      SCode.Encapsulated encflag;
      Absyn.Path fn;
      list<Absyn.Path> fns;
      FCore.Cache cache;
      InstanceHierarchy ih;
      SCode.Partial partialPrefix;
      Absyn.Info info;
      list<DAE.Function> resfns1,resfns2;
      SCode.Restriction rest;

    case (cache,_,ih,_,{}) then (cache,inEnv,ih,{});

    // Instantiate each function, add its FQ name to the type, needed when deoverloading
    case (cache,env,ih,_,(fn :: fns))
      equation
        // print("instOvl: " +& Absyn.pathString(fn) +& "\n");
        (cache,(c as SCode.CLASS(name=_,encapsulatedPrefix=_,restriction=rest)),cenv) =
           Lookup.lookupClass(cache, env, fn, true);
        true = SCode.isFunctionRestriction(rest);

        (cache,env,ih,resfns1) = implicitFunctionInstantiation2(inCache, cenv, inIH, DAE.NOMOD(), pre, c, {}, false);
        (cache,env,ih,resfns2) = instOverloadedFunctions(cache,env,ih,pre,fns);
      then (cache,env,ih,listAppend(resfns1,resfns2));

    // failure
    case (_,_,_,_,(fn :: _))
      equation
        Debug.fprint(Flags.FAILTRACE, "- Inst.instOverloaded_functions failed " +& Absyn.pathString(fn) +& "\n");
      then
        fail();
  end matchcontinue;
end instOverloadedFunctions;

protected function instExtDecl
"author: LS
  This function handles the external declaration. If there is an explicit
  call of the external function, the component references are looked up and
  inserted in the argument list, otherwise the input and output parameters
  are inserted in the argument list with their order. The return type is
  determined according to the specification; if there is a explicit call
  and a lhs, which must be an output parameter, the type of the function is
  that type. If no explicit call and only one output parameter exists, then
  this will be the return type of the function, otherwise the return type
  will be void."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input String inIdent;
  input SCode.ClassDef inClassDef;
  input Boolean inBoolean;
  input Prefix.Prefix inPrefix;
  input Absyn.Info info;
  output FCore.Cache outCache;
  output InnerOuter.InstHierarchy outIH;
  output DAE.ExternalDecl outExternalDecl;
algorithm
  (outCache,outIH,outExternalDecl) := matchcontinue (inCache,inEnv,inIH,inIdent,inClassDef,inBoolean,inPrefix,info)
    local
      String fname,lang,n;
      list<DAE.ExtArg> fargs;
      DAE.ExtArg rettype;
      Option<SCode.Annotation> ann;
      DAE.ExternalDecl daeextdecl;
      FCore.Graph env;
      SCode.ExternalDecl extdecl,orgextdecl;
      Boolean impl;
      list<SCode.Element> els;
      FCore.Cache cache;
      InstanceHierarchy ih;
      Prefix.Prefix pre;

    case (cache,env,ih,n,SCode.PARTS(elementLst=_,externalDecl = SOME(extdecl)),impl,pre,_) /* impl */
      equation
        InstUtil.isExtExplicitCall(extdecl);
        fname = InstUtil.instExtGetFname(extdecl, n);
        (cache,fargs) = InstUtil.instExtGetFargs(cache,env, extdecl, impl,pre,info);
        (cache,rettype) = InstUtil.instExtGetRettype(cache,env, extdecl, impl,pre,info);
        lang = InstUtil.instExtGetLang(extdecl);
        ann = InstUtil.instExtGetAnnotation(extdecl);
        daeextdecl = DAE.EXTERNALDECL(fname,fargs,rettype,lang,ann);
      then
        (cache,ih,daeextdecl);

    case (cache,env,ih,n,SCode.PARTS(elementLst = els,externalDecl = SOME(orgextdecl)),impl,pre,_)
      equation
        failure(InstUtil.isExtExplicitCall(orgextdecl));
        extdecl = InstUtil.instExtMakeExternaldecl(n, els, orgextdecl);
        (fname) = InstUtil.instExtGetFname(extdecl, n);
        (cache,fargs) = InstUtil.instExtGetFargs(cache,env, extdecl, impl,pre,info);
        (cache,rettype) = InstUtil.instExtGetRettype(cache,env, extdecl, impl,pre,info);
        lang = InstUtil.instExtGetLang(extdecl);
        ann = InstUtil.instExtGetAnnotation(orgextdecl);
        daeextdecl = DAE.EXTERNALDECL(fname,fargs,rettype,lang,ann);
      then
        (cache,ih,daeextdecl);
    else
      equation
        Debug.fprintln(Flags.FAILTRACE, "#-- Inst.instExtDecl failed");
      then
        fail();

  end matchcontinue;
end instExtDecl;

public function instRecordConstructorElt
"author: PA
  This function takes an Env and an Element and builds a input argument to
  a record constructor.
  E.g if the element is Real x; the resulting Var is \"input Real x;\""
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input InnerOuter.InstHierarchy inIH;
  input SCode.Element inElement;
  input DAE.Mod outerMod;
  input Boolean inImplicit;
  output FCore.Cache outCache;
  output InnerOuter.InstHierarchy outIH;
  output DAE.Var outVar;
algorithm
  (outCache,outIH,outVar) :=
  matchcontinue (inCache,inEnv,inIH,inElement,outerMod,inImplicit)
    local
      SCode.Element cl;
      FCore.Graph cenv,env,compenv;
      DAE.Mod mod_1;
      Absyn.ComponentRef owncref;
      DAE.Dimensions dimexp;
      DAE.Type tp_1;
      DAE.Binding bind;
      String id;
      SCode.Replaceable repl;
      SCode.Visibility vis;
      SCode.ConnectorType ct;
      Boolean impl;
      SCode.Attributes attr;
      list<Absyn.Subscript> dim;
      SCode.Parallelism prl;
      SCode.Variability var;
      Absyn.Direction dir;
      Absyn.Path t;
      SCode.Mod mod;
      SCode.Comment comment;
      SCode.Element elt;
      FCore.Cache cache;
      Absyn.InnerOuter io;
      SCode.Final finalPrefix;
      Absyn.Info info;
      InstanceHierarchy ih;
      Option<Absyn.ConstrainClass> cc;
      SCode.Prefixes prefixes;

    case (cache,env,ih,
          SCode.COMPONENT(name = id,
                          prefixes = prefixes as SCode.PREFIXES(
                            replaceablePrefix = _,
                            visibility = vis,

                            innerOuter = _
                          ),
                          attributes = (attr as
                          SCode.ATTR(arrayDims = dim, connectorType = ct,
                                     parallelism = prl,variability = var,direction = dir)),
                          typeSpec = Absyn.TPATH(t, _),modifications = mod,
                          comment = comment,
                          info = info),
          _,impl)
      equation
        // - Prefixes (constant, parameter, final, discrete, input, output, ...) of the remaining record components are removed.
        var = SCode.VAR();
        dir = Absyn.INPUT();
        attr = SCode.ATTR(dim,ct,prl,var,dir);

        (cache,cl,cenv) = Lookup.lookupClass(cache,env, t, true);
        (cache,mod_1) = Mod.elabMod(cache, env, ih, Prefix.NOPRE(), mod, impl, Mod.COMPONENT(id), info);
        mod_1 = Mod.merge(outerMod,mod_1,cenv,Prefix.NOPRE());
        owncref = Absyn.CREF_IDENT(id,{});
        (cache,dimexp) = InstUtil.elabArraydim(cache, env, owncref, t, dim, NONE(), false, NONE(), true, false, Prefix.NOPRE(), info, {});

        // cenv = FGraph.createVersionScope(env, cenv, SCode.elementName(cl), id, Prefix.NOPRE(), mod_1);
        (cache,compenv,ih,_,_,_,tp_1,_) = InstVar.instVar(cache, cenv, ih, UnitAbsyn.noStore, ClassInf.FUNCTION(Absyn.IDENT(""), false), mod_1, Prefix.NOPRE(),
          id, cl, attr, prefixes, dimexp, {}, {}, impl, comment, info, ConnectionGraph.EMPTY, Connect.emptySet, env);

        (cache,bind) = InstBinding.makeBinding(cache,env, attr, mod_1, tp_1, Prefix.NOPRE(), id, info);
      then
        (cache,ih,DAE.TYPES_VAR(id,DAE.ATTR(ct,prl,var,dir,Absyn.NOT_INNER_OUTER(),vis),tp_1,bind,NONE()));

    case (_,_,_,elt,_,_)
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.trace("- Inst.instRecordConstructorElt failed.,elt:");
        Debug.traceln(SCodeDump.unparseElementStr(elt,SCodeDump.defaultOptions));
      then
        fail();
  end matchcontinue;
end instRecordConstructorElt;

public function getRecordConstructorFunction
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input Absyn.Path inPath;
  output FCore.Cache outCache;
  output DAE.Function outFunc;
algorithm
  (outCache,outFunc)  := matchcontinue (inCache,inEnv,inPath)
    local
      Absyn.Path path;
      SCode.Element recordCl;
      FCore.Graph recordEnv;
      DAE.Function func;
      FCore.Cache cache;
      DAE.Type recType,fixedTy,funcTy;
      list<DAE.Var> vars, inputs, locals;
      list<DAE.FuncArg> fargs;
      DAE.EqualityConstraint eqCo;
      DAE.TypeSource src;
      String name, newName;

      case(_, _, _)
        equation
          path = Absyn.makeFullyQualified(inPath);
          func = FCore.getCachedInstFunc(inCache,path);
        then
          (inCache,func);

      case(_, _, _)
        equation
          (_,recordCl,recordEnv) = Lookup.lookupClass(inCache, inEnv, inPath, false);
          true = SCode.isRecord(recordCl);

          name = SCode.getElementName(recordCl);
          newName = FGraph.getInstanceOriginalName(recordEnv, name);
          recordCl = SCode.setClassName(newName, recordCl);

          (cache,_,_,_,_,_,recType,_,_,_) = Inst.instClass(inCache,recordEnv, InnerOuter.emptyInstHierarchy,
            UnitAbsynBuilder.emptyInstStore(), DAE.NOMOD(), Prefix.NOPRE(), recordCl,
            {}, true, InstTypes.INNER_CALL(), ConnectionGraph.EMPTY, Connect.emptySet);

          DAE.T_COMPLEX(ClassInf.RECORD(path), vars, eqCo, src) = recType;

          vars = Types.filterRecordComponents(vars, SCode.elementInfo(recordCl));
          (inputs,locals) = List.extractOnTrue(vars, Types.isModifiableTypesVar);
          inputs = List.map(inputs,Types.setVarDefaultInput);
          locals = List.map(locals,Types.setVarProtected);
          vars = listAppend(inputs,locals);

          path = Absyn.makeFullyQualified(path);
          fixedTy = DAE.T_COMPLEX(ClassInf.RECORD(path), vars, eqCo, src);
          fargs = Types.makeFargsList(inputs);
          funcTy = DAE.T_FUNCTION(fargs, fixedTy, DAE.FUNCTION_ATTRIBUTES_DEFAULT, {path});
          func = DAE.RECORD_CONSTRUCTOR(path,funcTy,DAE.emptyElementSource,DAE.VARIABLE());

          cache = InstUtil.addFunctionsToDAE(cache, {func}, SCode.NOT_PARTIAL());

          // add the instance record constructor too!
          path = Absyn.pathSetLastIdent(path, Absyn.makeIdentPathFromString(name));
          fixedTy = DAE.T_COMPLEX(ClassInf.RECORD(path), vars, eqCo, src);
          fargs = Types.makeFargsList(inputs);
          funcTy = DAE.T_FUNCTION(fargs, fixedTy, DAE.FUNCTION_ATTRIBUTES_DEFAULT, {path});
          func = DAE.RECORD_CONSTRUCTOR(path,funcTy,DAE.emptyElementSource,DAE.VARIABLE());

          cache = InstUtil.addFunctionsToDAE(cache, {func}, SCode.NOT_PARTIAL());

        then
          (cache,func);

      else
        equation
          true = Flags.isSet(Flags.FAILTRACE);
          Debug.fprint(Flags.FAILTRACE, "InstFunction.getRecordConstructorFunction failed for " +& Absyn.pathString(inPath) +& "\n");
        then
          fail();

  end matchcontinue;

end getRecordConstructorFunction;

public function addRecordConstructorFunction "Add record constructor whenever we instantiate a variable. Needed so we can cast to this constructor freely."
  input FCore.Cache inCache;
  input FCore.Graph inEnv;
  input DAE.Type inType;
  output FCore.Cache outCache;
algorithm
  outCache := matchcontinue (inCache,inEnv,inType)
    local
      list<DAE.Var> vars, inputs, locals;
      DAE.Type ty,recType,fixedTy,funcTy;
      DAE.EqualityConstraint eqCo;
      DAE.TypeSource src;
      FCore.Cache cache;
      Absyn.Path path;
      SCode.Element recordCl;
      FCore.Graph recordEnv;
      DAE.Function func;
      list<DAE.FuncArg> fargs;

    // try to instantiate class
    case (cache, _, DAE.T_COMPLEX(ClassInf.RECORD(path), _, _, _))
      equation
        path = Absyn.makeFullyQualified(path);
        (cache, _) = getRecordConstructorFunction(cache, inEnv, path);
      then
        cache;

    // if previous stuff didn't work, try to use the ty directly
    case (cache, _, DAE.T_COMPLEX(ClassInf.RECORD(path), vars, eqCo, src))
      equation
        path = Absyn.makeFullyQualified(path);

        (inputs,locals) = List.extractOnTrue(vars, Types.isModifiableTypesVar);
        inputs = List.map(inputs,Types.setVarDefaultInput);
        locals = List.map(locals,Types.setVarProtected);
        vars = listAppend(inputs,locals);

        fixedTy = DAE.T_COMPLEX(ClassInf.RECORD(path), vars, eqCo, src);
        fargs = Types.makeFargsList(inputs);
        funcTy = DAE.T_FUNCTION(fargs, fixedTy, DAE.FUNCTION_ATTRIBUTES_DEFAULT, {path});
        func = DAE.RECORD_CONSTRUCTOR(path,funcTy,DAE.emptyElementSource,DAE.VARIABLE());

        cache = InstUtil.addFunctionsToDAE(cache, {func}, SCode.NOT_PARTIAL());
      then
        (cache);

    else inCache;

  end matchcontinue;
end addRecordConstructorFunction;

protected function isElementImportantForFunction
  input SCode.Element elt;
  output Boolean b;
algorithm
  b := match elt
    case SCode.COMPONENT(prefixes=SCode.PREFIXES(visibility=SCode.PROTECTED()),
                         attributes=SCode.ATTR(direction=Absyn.BIDIR(),variability=SCode.VAR()))
      then false;
    else true;
  end match;
end isElementImportantForFunction;

protected function checkExtObjOutput
  input DAE.Type inType;
  input Absyn.Info info;
algorithm
  _ := match (inType,info)
    local
      Absyn.Path path;
      DAE.Type ty;
    case (DAE.T_FUNCTION(funcResultType=ty,source={path}),_)
      equation
        ((_,(_,_,true))) = Types.traverseType((ty,(path,info,true)),checkExtObjOutputWork);
      then ();
  end match;
end checkExtObjOutput;

protected function checkExtObjOutputWork
  input tuple<DAE.Type,tuple<Absyn.Path,Absyn.Info,Boolean>> inTpl;
  output tuple<DAE.Type,tuple<Absyn.Path,Absyn.Info,Boolean>> outTpl;
algorithm
  outTpl := match inTpl
    local
      Absyn.Path path1,path2;
      Absyn.Info info;
      String str1,str2;
      DAE.Type ty;
      Boolean b;
    case ((ty as DAE.T_COMPLEX(complexClassType=ClassInf.EXTERNAL_OBJ(path1)),(path2,info,true)))
      equation
        path1 = Absyn.joinPaths(path1,Absyn.IDENT("constructor"));
        str1 = Absyn.pathStringNoQual(path2);
        str2 = Absyn.pathStringNoQual(path1);
        b = Absyn.pathEqual(path1,path2);
        Error.assertionOrAddSourceMessage(b, Error.FUNCTION_RETURN_EXT_OBJ, {str1,str2}, info);
        outTpl = Util.if_(b,inTpl,(ty,(path2,info,false)));
      then outTpl;
    else inTpl;
  end match;
end checkExtObjOutputWork;

end InstFunction;
