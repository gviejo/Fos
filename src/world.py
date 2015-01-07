#!/usr/bin/env python
# -*- coding: utf-8 -*-


## @file world.py
# Gestion du "monde"  
#

import numpy as np

#-------------------------------------------
## Noeud du labyrinthe
#
# Peut être de nature variée : || , |_| , Y  , X , etc etc
class Noeud:
  #-----------------------------------------
  # Direction 0 is the North
  # Angles are positive counter-clockwise
  ##@var x
  # ( float )Abscisse du noeud
  ##@var y
  # ( float ) Ordonnée du noeud
  ##@var reward
  # ( float )  récompense associée au noeud
  ##@var length
  # 
  ##@var type
  # ( string )  type du noeud  : 'I', 'Y', 'i', 'T' 
  ##@var nbPxl
  # ( integer ) = CamDef (Définition de la "caméra" ).
  ##@var rot
  # ( float )  angle de rotation entre le noeud et le noeud standard de même type
  ##@var alleys
  # ( list of float ) liste des angles par rapport au nord des liens du noeud ( 3 liens pour un noeud de type "Y" par ex)
  ##@var goodChoice
  # ( boolean ) En allocentrique, cet attribut permet de savoir si le noeud étudié fait partit du plus court chemin entre le 
  # point de départ et la récompense. C'est le module de planification qui décide de la valeur du booléen.
  ##@var edges
  # ( list of  integer ) liste d'entiers correspondant aux identifiants des autres noeuds atteignables depuis le noeud étudié
  
  #-----------------------------------------
  ## Constructeur
  # @param rew recompense (reward) a associer au noeud
  # @param length ff
  # @param type 
  # @param rot rotation du noeud par rapport au noeud standard du même type (par défaut à 0°)
  # @param edges
  # @param camDef
  # @param goodChoice
  def __init__(self, rew, length, type, rot=0, edges=[],camDef=60, goodChoice=0):
    self.x = 0
    self.y = 0
    self.reward = rew
    self.length = length
    
    self.type = type

    self.nbPxl = camDef
    d = 360/camDef
    self.rot = rot/d
 
    # Positions of the alleys wrt the North, in pixels
    self.alleys = []
    self.goodChoice = goodChoice # The best alley choice to reach the platform, in allocentric coordinates

    if type == "I":
      self.alleys.append((0+self.rot)%self.nbPxl)
      self.alleys.append((180/d+self.rot)%self.nbPxl)
    elif type == "Y":
      self.alleys.append((54/d+self.rot)%self.nbPxl)
      self.alleys.append((180/d+self.rot)%self.nbPxl)
      self.alleys.append((306/d+self.rot)%self.nbPxl)
    elif type == "i":
      self.alleys.append((0+self.rot)%self.nbPxl)
    elif type == "T":
      self.alleys.append((90/d+self.rot)%self.nbPxl)
      self.alleys.append((180/d+self.rot)%self.nbPxl)
      self.alleys.append((270/d+self.rot)%self.nbPxl)
    elif type == "IY":
      self.alleys.append((54/d+self.rot)%self.nbPxl)
      self.alleys.append((180/d+self.rot)%self.nbPxl)
      #self.alleys.append((306/d+self.rot)%self.nbPxl)

    # edges: list of the names of the nodes reachable through each alley
    # the order is important for consistency
    if len(edges) != len(self.alleys):
      print "Node buiding pb: edges & alleys do not correspond."
    else:
      self.edges = edges

  #-----------------------------------------
  #def __del__(self):

  #-----------------------------------------
  ## Méthode pour placer le noeud ( modifie x et y )
  # @param x Abscisse du noeud
  # @param y Ordonne du noeud
  def setCoord(self,x,y):
    self.x, self.y = x, y

  #-----------------------------------------
  ## Simulate guidance of the mouse by the experimenter, returns the good choice
  # @param referentiel indique si on est en "allo" ou "ego" 
  # @param dir angle de la souris par rapport au nord EN PIXELS !
  def readGoodChoice(self, referentiel="allo", dir=0): # dir en pixels !
    if referentiel == "allo":
      return self.alleys[self.goodChoice]
    elif referentiel == "ego":
      #print "Classe Noeud: readGoodChoice:",referentiel,self.goodChoice,dir,self.nbPxl,(self.alleys[self.goodChoice]-dir)%self.nbPxl,self.alleys
      return (self.alleys[self.goodChoice]-dir)%self.nbPxl
    elif referentiel == "alloALaLouche":
      if self.alleys[self.goodChoice] == 0:
        return 0
      elif self.alleys[self.goodChoice]>0 and self.alleys[self.goodChoice]<self.nbPxl/2:
        return 1
      elif self.alleys[self.goodChoice] == self.nbPxl/2:
        return 2
      elif self.alleys[self.goodChoice]>self.nbPxl/2:
        return 3
    elif referentiel == "egoALaLouche":
      if (self.alleys[self.goodChoice]-dir)%self.nbPxl == 0:
        return 0
      elif (self.alleys[self.goodChoice]-dir)%self.nbPxl > 0 and self.alleys[self.goodChoice] < self.nbPxl/2:
        return 1
      elif (self.alleys[self.goodChoice]-dir)%self.nbPxl == self.nbPxl/2:
        return 2
      elif (self.alleys[self.goodChoice]-dir)%self.nbPxl > self.nbPxl/2:
        return 3
    else:
      print "Unkown referential:",referentiel    

  #-----------------------------------------
  # retourne le type de noeud
  def readType(self):
    type2vect = {'i':0,'I':1,'Y':2,'T':3,'IY':4}
    vect = np.zeros((len(type2vect)))
    vect[type2vect[self.type]] = 1
    return vect
    
  #-----------------------------------------
  ## Retourne la liste des angles correspondant aux allées disponibles pour la souris dans le noeud
  #  Les angles sont absolus en allocentrique, et relatif par rapport à l'angle de la souris en égocentrique
  def readCam(self, referentiel="allo", dir=0): # dir en pixels !
    cam = np.zeros((self.nbPxl))

    if referentiel == "allo":
      for i in self.alleys:
        cam[i]=1
    elif referentiel == "ego":
      for i in self.alleys:
        cam[(i-dir)%self.nbPxl]=1
        #cam[int(round((i-dir*self.nbPxl/360.)))%self.nbPxl]=1
        #print "- alley ", i, dir, (i-dir*self.nbPxl/360.), int(round((i-dir*self.nbPxl/360.))), (i-dir*self.nbPxl/360.)%self.nbPxl, cam
      #print "- inside readCam : ref : ", referentiel, "alleys : ",self.alleys, "dir : ", dir," camera : ", cam

    return cam

  #-----------------------------------------
  ## Lit la liste créée par readCam
  # Principalement pour faire du débug
  def displayCam(self, referentiel="allo", dir=0):
    cam = self.readCam(referentiel, dir)
    s = ""
    for i in cam:
      s = s + str(int(i))

    #print s

#-------------------------------------------
## Classe Labyrinthe
#
# Regroupe les noeuds créés et fait les liaisons.
class Maze:
  #-----------------------------------------
  ## @var type
  # ( string ) forme du labyrinthe
  ## @var Ns
  # ( dict ) dictionnaire des noeuds : {ID:data}
  ## @var DicoNoeuds
  # dictionnaire permettant le formattage du vecteur de sortie: {ID:indiceDsLeVecteur}
  # ça permet d'avoir toujours les mêmes noeuds dans le même ordre, car avec un iteritems, c'est pas gagné
  ## @var camDef
  # nombre d'actions possibles
  ## @var root
  # noeud racine du graphe, définition nécessaire pour les parcours de graphe
  
  #-----------------------------------------
  ## Constructeur
  # @param type forme du labyrinthe
  # @param camDef nombre d'actions possibles ( par défaut à 60 ) 
  def __init__(self, type, camDef=60):
    self.type = type       
    self.Ns = {}           
    self.DicoNoeuds = {}   

    self.camDef = camDef   
    self.root = "popo"    
    
    if type == "star":
      self.setStarStructure()
    elif type == "basicT":
      self.setBasicTStructure()
    elif type == "TYM":
      self.setTYMStructure()
    elif type == "2alleys":
      self.set2alleysStructure()
    else:
      print "Maze type unknown: ",type

  #-----------------------------------------
  ## Place les noeuds de manière récursive ( parcours DFS ) en partant du root.
  # @param rootID identifiant du noeud à placer
  # @param rootX abscisse du noeud
  # @param rootY ordonnée du noeud
  # @param setList liste des voisins du noeud sélectioné
  def setCoord(self, rootID, rootX, rootY, setList):
    self.Ns[rootID].setCoord(rootX,rootY)
    #print "Noeud ",rootID,"x",rootX,"y",rootY
    for i in range(len(self.Ns[rootID].edges)):
      cible = self.Ns[rootID].edges[i]
      if cible not in setList :
        setList.append(cible)
        offsetx = self.Ns[rootID].length * np.cos( (self.Ns[rootID].alleys[i] +self.Ns[rootID].nbPxl/4.) / float(self.Ns[rootID].nbPxl) *2*np.pi)
        offsety = self.Ns[rootID].length * np.sin( (self.Ns[rootID].alleys[i] +self.Ns[rootID].nbPxl/4.) / float(self.Ns[rootID].nbPxl) *2*np.pi)
        self.setCoord(cible, rootX+offsetx, rootY+offsety, setList)
 
  #-----------------------------------------
  ## Donne l'ID du noeud ayant une récompense de 1
  def findRew(self):
    for ID,n in self.Ns.iteritems():
      if self.Ns[ID].reward == 1 :
        return ID

  #-----------------------------------------
  ## Construit le labyrinthe en "T"
  def setBasicTStructure(self):
    self.root = "0"

    self.Ns["0"]=Noeud(0,1,"i",0,["1"],self.camDef)
    self.Ns["1"]=Noeud(0,1,"T",0,["2","0","3"],self.camDef)
    self.Ns["2"]=Noeud(0,1,"i",-90,["1"],self.camDef)
    self.Ns["3"]=Noeud(1,1,"i",90,["1"],self.camDef)

    self.makeDicoNoeuds()

    self.setCoord(self.root,0,0,[self.root])
 
    self.nbTransitionsMax = 10
 
  #-----------------------------------------
  ## Construit le starmaze
  def setStarStructure(self):
    self.root = "1a"

    self.Ns["I"]=Noeud(0,0.235,"Y",0,["10","1b","2"],goodChoice=0)

    self.Ns["1a"]=Noeud(0,0.235,"i",0,["1b"],goodChoice=0)
    self.Ns["1b"]=Noeud(0,0.235,"I",0,["I","1a"],goodChoice=0)

    self.Ns["2"]=Noeud(0,0.235,"I",7*18,["I","V"],goodChoice=0)
    
    self.Ns["V"]=Noeud(0,0.235,"Y",4*18,["2","3b","4"],goodChoice=2) # quel guidage est fait dans ce cas ?

    self.Ns["3a"]=Noeud(0,0.235,"i",4*18,["3b"],goodChoice=0)
    self.Ns["3b"]=Noeud(0,0.235,"I",4*18,["V","3a"],goodChoice=0)

    self.Ns["4"]=Noeud(0,0.235,"I",11*18,["V","IV"],goodChoice=1)

    self.Ns["IV"]=Noeud(0,0.235,"Y",8*18,["4","5b","6"],goodChoice=2)

    self.Ns["5a"]=Noeud(0,0.235,"i",8*18,["5b"],goodChoice=0)
    self.Ns["5b"]=Noeud(0,0.235,"I",8*18,["IV","5a"],goodChoice=0)

    self.Ns["6"]=Noeud(0,0.235,"I",15*18,["IV","III"],goodChoice=1)

    self.Ns["III"]=Noeud(0,0.235,"Y",12*18,["6","7b","8"],goodChoice=1)

    self.Ns["7a"]=Noeud(1,0.235,"i",12*18,["7b"],goodChoice=0)
    self.Ns["7b"]=Noeud(0,0.235,"I",12*18,["III","7a"],goodChoice=1)

    self.Ns["8"]=Noeud(0,0.235,"I",19*18,["III","II"],goodChoice=0)

    self.Ns["II"]=Noeud(0,0.235,"Y",16*18,["8","9b","10"],goodChoice=0)

    self.Ns["9a"]=Noeud(0,0.235,"i",16*18,["9b"],goodChoice=0)
    self.Ns["9b"]=Noeud(0,0.235,"I",16*18,["II","9a"],goodChoice=0)
    
    self.Ns["10"]=Noeud(0,0.235,"I",3*18,["II","I"],goodChoice=0)

    self.makeDicoNoeuds()

    self.setCoord(self.root,0,0,[self.root])
    
    self.nbTransitionsMax = 57

  #-----------------------------------------
  ## Construit le TYM
  def setTYMStructure(self):
    self.root = "1a"

    self.Ns["I"]=Noeud(0,0.235,"Y",0,["10","1b","2"],goodChoice=0)

    self.Ns["1a"]=Noeud(0,0.235,"i",0,["1b"],goodChoice=0)
    self.Ns["1b"]=Noeud(0,0.235,"I",0,["I","1a"],goodChoice=0)

    self.Ns["2"]=Noeud(0,0.235,"I",7*18,["I","V"],goodChoice=0)
    
    self.Ns["V"]=Noeud(0,0.235,"Y",4*18,["2","3b","4"],goodChoice=0) # quel guidage est fait dans ce cas ?

    self.Ns["3a"]=Noeud(0,0.235,"i",4*18,["3b"],goodChoice=0)
    self.Ns["3b"]=Noeud(0,0.235,"I",4*18,["V","3a"],goodChoice=0)

    self.Ns["4"]=Noeud(0,0.235,"I",11*18,["V","IV"],goodChoice=0)

    self.Ns["IV"]=Noeud(0,0.235,"i",11*18,["4"],goodChoice=0)

    self.Ns["III"]=Noeud(1,0.235,"i",9*18,["8"],goodChoice=0)

    self.Ns["8"]=Noeud(0,0.235,"I",19*18,["III","II"],goodChoice=0)

    self.Ns["II"]=Noeud(0,0.235,"Y",16*18,["8","9b","10"],goodChoice=0)

    self.Ns["9a"]=Noeud(0,0.235,"i",16*18,["9b"],goodChoice=0)
    self.Ns["9b"]=Noeud(0,0.235,"I",16*18,["II","9a"],goodChoice=0)

    self.Ns["10"]=Noeud(0,0.235,"I",3*18,["II","I"],goodChoice=0)

    self.makeDicoNoeuds()

    self.setCoord(self.root,0,0,[self.root])
    
    self.nbTransitionsMax = 38

  #-----------------------------------------

  ## Construit la version 2 couloirs pour les contrôles nage Fos
  def set2alleysStructure(self):
    self.root = "1a"

    self.Ns["I"]=Noeud(0,0.235,"IY",0,["10","1b"],goodChoice=0)

    self.Ns["1a"]=Noeud(0,0.235,"i",0,["1b"],goodChoice=0)
    self.Ns["1b"]=Noeud(0,0.235,"I",0,["I","1a"],goodChoice=0)

    self.Ns["II"]=Noeud(0,0.235,"i",13*18,["10"],goodChoice=0)

    self.Ns["10"]=Noeud(0,0.235,"I",3*18,["II","I"],goodChoice=0)

    self.makeDicoNoeuds()

    self.setCoord(self.root,0,0,[self.root])
    
    self.nbTransitionsMax = 38


  #-----------------------------------------

  ## Ancienne version du starmaze
  ## Ancienne version du starmaze
  def setOldStarStructure(self):
    self.root = "1a"

    self.Ns["I"]=Noeud(0,0.235,"Y",0,["10b","1b","2a"])

    self.Ns["1a"]=Noeud(0,0.235,"i",0,["1b"])
    self.Ns["1b"]=Noeud(0,0.235,"I",0,["I","1a"])

    self.Ns["2a"]=Noeud(0,0.235,"I",7*18,["I","2b"])
    self.Ns["2b"]=Noeud(0,0.235,"I",7*18,["2a","V"])
    
    self.Ns["V"]=Noeud(0,0.235,"Y",4*18,["2b","3b","4a"])

    self.Ns["3a"]=Noeud(0,0.235,"i",4*18,["3b"])
    self.Ns["3b"]=Noeud(0,0.235,"I",4*18,["V","3a"])

    self.Ns["4a"]=Noeud(0,0.235,"I",11*18,["V","4b"])
    self.Ns["4b"]=Noeud(0,0.235,"I",11*18,["4a","IV"])

    self.Ns["IV"]=Noeud(0,0.235,"Y",8*18,["4b","5b","6a"])

    self.Ns["5a"]=Noeud(0,0.235,"i",8*18,["5b"])
    self.Ns["5b"]=Noeud(0,0.235,"I",8*18,["IV","5a"])

    self.Ns["6a"]=Noeud(0,0.235,"I",15*18,["IV","6b"])
    self.Ns["6b"]=Noeud(0,0.235,"I",15*18,["6a","III"])

    self.Ns["III"]=Noeud(0,0.235,"Y",12*18,["6b","7b","8a"])

    self.Ns["7a"]=Noeud(1,0.235,"i",12*18,["7b"])
    self.Ns["7b"]=Noeud(0,0.235,"I",12*18,["III","7a"])

    self.Ns["8a"]=Noeud(0,0.235,"I",19*18,["III","8b"])
    self.Ns["8b"]=Noeud(0,0.235,"I",19*18,["8a","II"])

    self.Ns["II"]=Noeud(0,0.235,"Y",16*18,["8b","9b","10a"])

    self.Ns["9a"]=Noeud(0,0.235,"i",16*18,["9b"])
    self.Ns["9b"]=Noeud(0,0.235,"I",16*18,["II","9a"])

    self.Ns["10a"]=Noeud(0,0.235,"I",3*18,["II","10b"])
    self.Ns["10b"]=Noeud(0,0.235,"I",3*18,["10a","I"])

    self.setCoord(self.root,0,0,[self.root])

    self.makeDicoNoeuds()
    
    self.nbTransitionsMax = 57

  #-----------------------------------------
  ## Construit le dictionnaire des identifiants des noeuds du labyrinthe
  def makeDicoNoeuds(self):
    
    i = 0
    for ID,n in self.Ns.iteritems():
      self.DicoNoeuds[ID]=i
      i+=1

    #print self.DicoNoeuds

#-------------------------------------------
## En plus de contenir la description de l'environnement (le graphe de Maze),
 # il y a également les infos sur la position de la souris, etc.
class World:
  #-----------------------------------------
  ## @var M
  # ( Maze ) labyrinthe
  ## @var mousePos
  # ( integer) ID du noeud sur lequel est la souris présentement
  ## @var mouseDir
  # ( integer) angle de la souris par rapport au nord
  ## @var t
  # ( integer) variable de temps
  ## @var f
  # ( file ) variable "file" pour les logs
  
  
  #-----------------------------------------
  ## Constructeur
  # @param type forme du labyrinthe
  # @param camDef nombre d'actions possibles ( par défaut à 60 ) 
  def __init__(self, type, camDef=60):
    self.M = Maze(type, camDef)

    self.mousePos = self.M.root
    self.mouseDir = 0

    self.t = 0

    # self.f=open('log/World','w')    # log file where the internal state will be stored if logAll function is used
    # self.log()

  #-----------------------------------------
  ## Destructeur
  # def __del__(self):
  #   self.f.close()

  #-----------------------------------------
  ## Place la souris en position initiale
  # @param pos Identifiant ( integer ) du noeud où la souris doit partir ( root par défaut)
  # @param dir direction de la souris par rapport au nord ( 0 par defaut)
  def startingPos(self,pos=-1,dir=0):
    if pos == -1:
      self.mousePos = self.M.root
    else :
      self.mousePos = pos
    self.mouseDir = dir

  #-----------------------------------------
  ## to simulate guidance of the mouse by the experimenter, returns the good choice
  # @param referentiel referentiel de la souris ( 'allo' pour allocentrique, 'ego' pour égocentrique)
  def readGoodChoice(self, referentiel="allo"):
    #print "W : readgoodchoice : ", self.mousePos, referentiel, self.mouseDir, "choix",self.M.Ns[self.mousePos].readGoodChoice(referentiel,self.mouseDir)
    return self.M.Ns[self.mousePos].readGoodChoice(referentiel,self.mouseDir)

  #-----------------------------------------
  ## renvoie un tableau de la taille du nb de noeuds, avec des 0 partout sauf
  # un 1 correspondant à l'emplacement courant de la souris
  # Utilisé par l'Acteur-Critique de type PRTR
  # @param pos Identifiant ( integer ) du noeud où la souris doit partir ( root par défaut)
  def readPlaceCells(self,pos=-1):
    PC = np.zeros((len(self.M.DicoNoeuds)))
    if pos ==-1:
      PC[self.M.DicoNoeuds[self.mousePos]] = 1
    else:
      PC[self.M.DicoNoeuds[pos]] = 1
    #print PC
    return PC

  #-----------------------------------------
  ## renvoie l'identifiant de la PC courante
  def readPlaceCellsID(self):
    #return self.M.DicoNoeuds[self.mousePos]
    return self.mousePos

  #-----------------------------------------
  ## renvoie la longueur du Noeud courant
  def readNodeLength(self):
    return self.M.Ns[self.mousePos].length

  #-----------------------------------------
  ##Retourne le type du noeud courant, utilisé pour les apprentissages sur une description grossière des états.
  def readType(self):
    return self.M.Ns[self.mousePos].readType()

  #----------------------------------------- 
  # retourne un vecteur de 4 composante indiquant s'il est possible
  # d'aller tout droit, de tourner à gauche, de faire demi-tour et de
  # tourner à droite
  def readPathwaysALaLouche(self,referentiel='ego'):
    pathways = np.zeros(4)
    if referentiel == 'ego':
      if self.mouseDir%self.M.camDef in self.M.Ns[self.mousePos].alleys:
        pathways[0]=1
      for i in range(1,self.M.camDef/2):
        if ((i+self.mouseDir)%self.M.camDef) in self.M.Ns[self.mousePos].alleys:
          pathways[1]=1
      if (self.M.camDef/2+self.mouseDir)%self.M.camDef in self.M.Ns[self.mousePos].alleys:
        pathways[2]=1
      for i in range(self.M.camDef-1,self.M.camDef/2,-1):
        if ((i+self.mouseDir)%self.M.camDef) in self.M.Ns[self.mousePos].alleys:
          pathways[3]=1
    elif referentiel == 'allo':
      if 0 in self.M.Ns[self.mousePos].alleys:
        pathways[0]=1
      for i in range(1,self.M.camDef/2):
        if i in self.M.Ns[self.mousePos].alleys:
          pathways[1]=1
      if self.M.camDef/2 in self.M.Ns[self.mousePos].alleys:
        pathways[2]=1
      for i in range(self.M.camDef-1,self.M.camDef/2,-1):
        if i in self.M.Ns[self.mousePos].alleys:
          pathways[3]=1

    return pathways

  #-----------------------------------------
  ##Retourne la liste des angles correspondant aux allées disponibles pour la souris dans le noeud.
  #Les angles sont absolus en allocentrique, et relatif par rapport à l'angle de la souris en égocentrique.
  def readCam(self, referentiel="allo"):
    return self.M.Ns[self.mousePos].readCam(referentiel,self.mouseDir)

  #-----------------------------------------
  ##Retourne la valeur de récompense associée à la position courante
  def readRew(self):
    return self.M.Ns[self.mousePos].reward

  #-----------------------------------------
  ## Move the mouse (if possible) according to very general orders : forward, backward, left or right
  # @param dir angle de la souris par rapport au nord
  def moveALaLouche(self,act,referentiel="ego"):
    moved = False
    directionsToTry = []
    if referentiel == "ego":
      if act == 0: # FORWARD
        directionsToTry.append(self.mouseDir%self.M.camDef)
        #print "Forward"
      elif act == 1: # LEFT
        for i in range(1,self.M.camDef/2):
          directionsToTry.append((i+self.mouseDir)%self.M.camDef)
        #print "Left"
      elif act == 2: # BACKWARD
        directionsToTry.append((self.M.camDef/2+self.mouseDir)%self.M.camDef)
        #print "Backward"
      elif act == 3: # RIGHT
        for i in range(self.M.camDef-1,self.M.camDef/2,-1):
          directionsToTry.append((i+self.mouseDir)%self.M.camDef)
        #print "W.moveALaLouche():Right",directionsToTry
    elif referentiel == "allo":
      if act == 0: # FORWARD
        directionsToTry.append(0)
        #print "Forward"
      elif act == 1: # LEFT
        for i in range(1,self.M.camDef/2):
          directionsToTry.append(i)
        #print "Left"
      elif act == 2: # BACKWARD
        directionsToTry.append(self.M.camDef/2)
        #print "Backward"
      elif act == 3: # RIGHT
        for i in range(self.M.camDef-1,self.M.camDef/2,-1):
          directionsToTry.append(i)
        #print "Right"

    for d in directionsToTry:
      if not moved and (d in self.M.Ns[self.mousePos].alleys):
        #DEBUG
        #if act==3:
        #  print "W.moveALaLouche():Right: dir=",self.mouseDir,d
        #  print "W.moveALaLouche():Right: pos=",self.mousePos
        #  print "W.moveALaLouche():Right: node=",self.M.Ns[self.mousePos].alleys,self.M.Ns[self.mousePos].edges
        #END DEBUG
        self.mouseDir = d
        self.mousePos = self.M.Ns[self.mousePos].edges[self.M.Ns[self.mousePos].alleys.index(d)]
        #print "Mouse moved in direction ", self.mouseDir, "arrived in ", self.mousePos
        moved = True

      
    if not moved:
      print "Mouse bumped in a wall !"
      print referentiel
      print "act=",act
      print "pos=",self.mousePos
      print "dir=",self.mouseDir
      exit()

    self.t += 1
    # self.log()
    return moved

  #-----------------------------------------
  ## Move the mouse (if possible) according to the direction dir (in pixels, in egocentric ref. frame)
  # @param dir angle de la souris par rapport au nord
  # @param referentiel referentiel de la souris ( 'allo' pour allocentrique, 'ego' pour égocentrique)
  def move(self,dir,referentiel="allo"):
    moved = False
    for i in range(len(self.M.Ns[self.mousePos].alleys)):
      #print i,dir+self.mouseDir, len(self.M.Ns[self.mousePos].alleys)
      if not moved :
        if referentiel == 'ego':
          if ((dir+self.mouseDir)%self.M.camDef == self.M.Ns[self.mousePos].alleys[i]):
            self.mouseDir = (dir+self.mouseDir)%self.M.camDef
            #print "edges ",self.M.Ns[self.mousePos].edges, self.M.Ns[self.mousePos].edges[i]
            self.mousePos = self.M.Ns[self.mousePos].edges[i]
            moved = True
            #print "Mouse moved to ", self.mousePos, "in direction ", self.mouseDir
        elif referentiel == 'allo':
          if (dir == self.M.Ns[self.mousePos].alleys[i]):
            self.mouseDir = dir
            #print "edges ",self.M.Ns[self.mousePos].edges, self.M.Ns[self.mousePos].edges[i]
            self.mousePos = self.M.Ns[self.mousePos].edges[i]
            moved = True
            #print "Mouse moved to ", self.mousePos, "in direction ", self.mouseDir          

    if not moved:
      print "Mouse bumped in a wall !"

    self.t += 1
    # self.log()
    return moved

  #-----------------------------------------
  ## Ecrit dans les logs selon la macro suivante : temps	mousePos	mouseDir
  def log(self):
    self.f.writelines(str(self.t)+" "+str(self.mousePos)+" "+str(self.mouseDir)+"\n")

#---------------------------
## Main du fichier
# Instancie un labyrinthe et place la souris en position initiale.
# Usage essentiellement réservé à du déboggage.
def main():
  # Test Noeud
  #-----------

  #N = Noeud(0,0.235,"Y",0)
  #N.displayCam()
  #print N.readType()

  # Test Maze
  #-----------
  #M = Maze("star")
  #M = Maze("basicT")
  #print M.Ns["I"].alleys
  #print "0    -    |    -    |    -    |    -    |    -    |    -    |"
  #M.Ns["I"].displayCam()
  #M.Ns["I"].displayCam("ego",0)
  #M.Ns["I"].displayCam("ego",90)

  # Test Maze
  #-----------
  #M = Maze("TYM")
  #M = Maze("basicT")
  #print M.Ns["I"].alleys
  #print "0    -    |    -    |    -    |    -    |    -    |    -    |"
  #M.Ns["I"].displayCam()
  #M.Ns["I"].displayCam("ego",0)
  #M.Ns["I"].displayCam("ego",90)

  # Test World
  #-----------

  W = World("TYM")
  W.moveALaLouche(0) # fwd
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(0) # fwd
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(0) # bump
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(1) # left
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(0) # fwd
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(3) # right
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

  W.moveALaLouche(0) # fwd
  print W.readPathwaysALaLouche()
  print W.readGoodChoice('alloALaLouche')
  print W.readGoodChoice('egoALaLouche')

#---------------------------

if __name__ == '__main__':
  # Import Psyco if available
  try:
    import psyco
    psyco.log()
    psyco.profile()
    psyco.full()
  except ImportError:
    print 'Psyco not available.'
  main()
