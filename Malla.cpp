// =================================================================================================
// ======================================================================================= Malla.cpp
//
// Malla d-dimensional de polihedros arbitrarios
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================== MÓDULO GLOBAL

module;

#include <gmsh.h>

// ===================================================================== DECLARACIÓN DE LA PARTICIÓN
// =================================================================================================

export module VF:Malla;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;

import std;

namespace VF
{
// ============================================================================== VARIABLES EN LÍNEA
// =========================================================================================== MaxID

inline constexpr std::size_t MaxID = std::numeric_limits<std::size_t>::max();

// =========================================================================== DECLARACIÓN DE CLASES
// ========================================================================================== TPunto

export
template<std::size_t d>
class TPunto : public TVector<d>
{
public:
  std::size_t const ID;                          // Identificador

// --------------------------------------------------------------------------------------- Funciones

public:
  TPunto() = delete;

  TPunto(std::size_t const ID_, std::convertible_to<double> auto const... Coord) :
    TVector<d>{Coord...}, ID(ID_) {}
};

// =================================================================================================
// =========================================================================================== TCara

export
template<std::size_t d>
class TCara
{
private:
  TCelda<d> const * const CeldaPPtr;              // Puntero a la celda propietaria
  TCelda<d> const *CeldaNPtr = nullptr;           // Puntero a la celda vecina
  std::vector<TPunto<d> const *> const PtoPtrVec; // Punteros a puntos en orden antihorario

public:
  TVector<d> const Sf,                            // Área
                   Cf;                            // Centroide
  std::size_t CCID = MaxID;                       // Identificador del contorno al que pertenece

// --------------------------------------------------------------------------------------- Funciones

private:
  TVector<d>
  Sf_() const;

  TVector<d>
  Cf_() const;

public:
  TCara() = delete;

  TCara(TCelda<d> const * const CeldaPPtr_, std::vector<TPunto<d> const *> &&PtoPtrVec_) :
    CeldaPPtr(CeldaPPtr_), PtoPtrVec(std::move(PtoPtrVec_)), Sf(Sf_()), Cf(Cf_()) {}

  TCelda<d> const
  &CeldaP() const
    { return *CeldaPPtr; }

  TCelda<d> const
  &CeldaN() const
    { return *CeldaNPtr; }

  std::size_t
  NPunto() const
    { return PtoPtrVec.size(); }

  TPunto<d> const
  &Punto(std::size_t const i) const
    { return *PtoPtrVec[i]; }

  bool
  EsCC() const
    { return CCID != MaxID; }

  TVector<d>
  nf() const
    { return Sf / mag(Sf); }

  TVector<d> const
  L() const
    { return (((Cf - CeldaP().C) & Sf) / (Sf & Sf)) * Sf; }

  TVector<d>
  δ() const
    { return CeldaN().C - CeldaP().C; }

  template<std::size_t r>
  TTensor<d, r>
  Interpola(TTensor<d, r> const &φP, TTensor<d, r> const &φN) const
    { return lerp(φP, φN, ((Cf - CeldaP().C) & Sf) / (δ() & Sf)); }

// ------------------------------------------------------------------------------------------ Amigos

  friend class TMalla<d>;
};

// =================================================================================================
// ========================================================================================== TCelda

export
template<std::size_t d>
class TCelda
{
public:
  enum ETipo { TRIAN, CUADR, TETRA, HEXA, CUÑA, PIRAM };

// ------------------------------------------------------------------------------------------- Datos

public:
  ETipo const Tipo;                               // Tipo de celda
  std::size_t const ID;                           // Identificador de la celda
  std::size_t const GrpID;                        // Identificador del grupo al que pertenece

private:
  std::vector<TPunto<d> const *> const PtoPtrVec; // Punteros a puntos en orden MSH/VTK
  std::vector<TCara<d>> CaraVec;                  // Caras

public:
  double const V;                                 // Volumen
  TVector<d> const C;                             // Centroide

  static constexpr std::size_t MaxNCara = 2u * d; // Número máximo de caras

// --------------------------------------------------------------------------------------- Funciones

private:
  template<ETipo>
  std::vector<TCara<d>>
  CaraVec_() const = delete;

  std::vector<TCara<d>>
  CaraVec_() const;

  double
  V_() const;

  TVector<d>
  C_() const;

public:
  TCelda() = delete;

  TCelda(ETipo const Tipo_, std::size_t ID_, std::size_t const GrpID_,
         std::vector<TPunto<d> const *> &&PtoPtrVec_) :
    Tipo(Tipo_), ID(ID_), GrpID(GrpID_), PtoPtrVec(std::move(PtoPtrVec_)), CaraVec(CaraVec_()),
    V(V_()), C(C_()) {}

  std::size_t
  NPunto() const
    { return PtoPtrVec.size(); }

  TPunto<d> const
  &Punto(std::size_t const i) const
    { return *PtoPtrVec[i]; }

  std::size_t
  NCara() const
    { return CaraVec.size(); }

  TCara<d> const
  &Cara(std::size_t const i) const
    { return CaraVec[i]; }

// ------------------------------------------------------------------------------------------ Amigos

  auto
  friend begin(TCelda const &Celda)
    { return std::begin(Celda.CaraVec); }

  auto
  friend begin(TCelda &Celda)
    { return std::begin(Celda.CaraVec); }

  auto
  friend end(TCelda const &Celda)
    { return std::end(Celda.CaraVec); }

  auto
  friend end(TCelda &Celda)
    { return std::end(Celda.CaraVec); }
};

// =================================================================================================
// ========================================================================================== TMalla

export
template<std::size_t d>
class TMalla
{
private:
  using TIDVec = std::array<std::size_t, TCelda<d>::MaxNCara>;

  class THelper
  {
  private:
    std::vector<std::size_t> TagIDVec;
    std::vector<std::vector<TCelda<d> *>> CeldaPtrVecVec;

  public:
    void
    DefNTag(std::size_t const NTag)
      { TagIDVec.resize(NTag); }

    void
    DefTag(std::size_t const ID, std::size_t const Tag)
      { TagIDVec[Tag - 1u] = ID; }

    std::size_t
    ID(std::size_t const Tag) const
      { return TagIDVec[Tag - 1u]; }

    void
    DefNPunto(std::size_t const NPunto)
      { CeldaPtrVecVec.resize(NPunto); }

    TPunto<d> const
    *PuntoPtr(std::size_t const Tag) const
      { return &TMalla<d>::Punto(ID(Tag)); }

    void
    DefCeldaPtr(std::size_t const Tag, TCelda<d> *CeldaPtr)
      { CeldaPtrVecVec[ID(Tag)].push_back(CeldaPtr); }

    std::vector<TCelda<d> *> const
    &CeldaPtrVec(std::size_t const ID) const
      { return CeldaPtrVecVec[ID]; }

    template<std::convertible_to<std::size_t>... TRest>
    bool
    EsCeldaN(TCelda<d> const *CeldaPtr, std::size_t const, TRest const... ID) const
      { return (... && std::ranges::contains(CeldaPtrVec(ID), CeldaPtr)); }
  };

  template<std::size_t = d, std::size_t = 1u>
  class TC : public TExprBase<TC<d, 1u>>
  {
  public:
    TVector<d> const
    &Eval(TCelda<d> const &Celda) const
      { return Celda.C; }

    TVector<d> const
    &Eval(TCara<d> const &Cara) const
      { return Cara.Cf; }

    template<bool>
    TVector<d> const
    &Coef(TCelda<d> const &Celda) const
      { return Eval(Celda); }
  };

  template<std::size_t = d, std::size_t = 0u>
  class TV : public TExprBase<TV<d, 0u>>
  {
  public:
    double
    Eval(TCelda<d> const &Celda) const
      { return Celda.V; }

    template<bool>
    double
    Coef(TCelda<d> const &Celda) const
      { return Eval(Celda); }
  };

// ------------------------------------------------------------------------------------------- Datos

private:
  static inline std::vector<TPunto<d>> PuntoVec;    // Puntos
  static inline std::vector<TCelda<d>> CeldaVec;    // Celdas
  static inline std::vector<TIDVec> IDVecVec;       // ID de las celdas vecinas
  static inline std::vector<std::string> GrpStrVec, // Grupos de celdas
                                         CCStrVec;  // Grupos de caras del contorno

// --------------------------------------------------------------------------------------- Funciones

private:
  void
  static ReadNodes(THelper &);

  template<int>
  void
  static ReadPhysGrp(std::vector<std::string> &, std::map<int, std::size_t> &);

  template<std::size_t... i>
  void
  static ReadElements(std::index_sequence<i...>, THelper &);

  template<std::size_t... i>
  void
  static ReadBoundary(std::index_sequence<i...>, THelper const &);

  template<std::size_t... i>
  void
  static ReadPeriodics(std::index_sequence<i...>, THelper const &);

  template<std::size_t... i>
  TPunto<d>
  static &DefPunto(std::index_sequence<i...>, std::span<double const> const Coord)
    { return PuntoVec.emplace_back(PuntoVec.size(), Coord[i]...); }

  TPunto<d>
  static &DefPunto(std::span<double const> const Coord)
    { return DefPunto(std::make_index_sequence<d>{}, Coord); }

  TCelda<d>
  static &DefCelda(TCelda<d>::ETipo const Tipo, std::size_t const GrpID,
                   std::vector<TPunto<d> const *> &&PtoPtrVec)
    { return CeldaVec.emplace_back(Tipo, CeldaVec.size(), GrpID, std::move(PtoPtrVec)); }

  std::size_t
  static ID(std::vector<std::string> const &, std::string_view const);

public:
  TMalla() = delete;

  void
  static Read(std::string const &);

  void
  static Write(std::ostream &);

  std::size_t
  static NPunto()
    { return PuntoVec.size(); }

  TPunto<d> const
  static &Punto(std::size_t const i)
    { return PuntoVec[i]; }

  std::size_t
  static NCelda()
    { return CeldaVec.size(); }

  TCelda<d> const
  static &Celda(std::size_t const i)
    { return CeldaVec[i]; }

  TIDVec const
  static &IDVec(std::size_t const i)
    { return IDVecVec[i]; }

  std::size_t
  static NCC()
    { return CCStrVec.size(); }

  std::size_t
  static CCID(std::string_view const CCStr)
    { return ID(CCStrVec, CCStr); }

  std::size_t
  static GrpID(std::string_view const GrpStr)
    { return ID(GrpStrVec, GrpStr); }

  TC<>
  static C()
    { return {}; }

  TV<>
  static V()
    { return {}; }

  auto
  static Begin()
    { return std::cbegin(CeldaVec); }

  auto
  static End()
    { return std::cend(CeldaVec); }
};

// =========================================================================================== ALIAS
// ======================================================================================== TMalla2D

export
using TMalla2D = TMalla<2u>;

// =================================================================================================
// ======================================================================================== TMalla3D

export
using TMalla3D = TMalla<3u>;

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// =========================================================================================== TCara

template<>
TVector<2u>
inline TCara<2u>::Sf_() const
{
TVector<2u> const r = Punto(1u) - Punto(0u);

return {r[1u], -r[0u]};
}

// =================================================================================================

template<>
TVector<3u>
inline TCara<3u>::Sf_() const
{
TVector<3u> Sf = {};

for (std::size_t i = 0u; i < NPunto(); ++i)
  Sf += Punto(i) ^ Punto((i + 1u) % NPunto());
return 0.5 * Sf;
}

// =================================================================================================

template<>
TVector<2u>
inline TCara<2u>::Cf_() const
  { return 0.5 * (Punto(0u) + Punto(1u)); }

// =================================================================================================

template<>
TVector<3u>
inline TCara<3u>::Cf_() const
{
TPunto<3u> const &Pto = Punto(0u);
TVector<3u> Cf = {};

for (std::size_t i = 2u; i < NPunto(); ++i)
  Cf += (Punto(i - 1u) + Punto(i)) * mag((Punto(i - 1u) - Pto) ^ (Punto(i) - Pto));
Cf /= 6.0 * mag(Sf);
return Cf + Pto / 3.0;
}

// =================================================================================================
// ========================================================================================== TCelda

template<>                                                      /*    2                           */
template<>                                                      /*    |`\                         */
std::vector<TCara<2u>>                                          /*    |  `\                       */
inline TCelda<2u>::CaraVec_<TCelda<2u>::TRIAN>() const          /*    |    `\                     */
{                                                               /*    |      `\                   */
TVector const u = *PtoPtrVec[1u] - *PtoPtrVec[0u],              /*    |        `\                 */
              v = *PtoPtrVec[2u] - *PtoPtrVec[0u];              /*    0----------1                */

if (u[0u] * v[1u] > u[1u] * v[0u])
  return {{this, {PtoPtrVec[0u], PtoPtrVec[1u]}},
          {this, {PtoPtrVec[1u], PtoPtrVec[2u]}},
          {this, {PtoPtrVec[2u], PtoPtrVec[0u]}}};
else
  return {{this, {PtoPtrVec[0u], PtoPtrVec[2u]}},
          {this, {PtoPtrVec[2u], PtoPtrVec[1u]}},
          {this, {PtoPtrVec[1u], PtoPtrVec[0u]}}};
}

// =================================================================================================

template<>                                                      /*    3-----------2               */
template<>                                                      /*    |           |               */
std::vector<TCara<2u>>                                          /*    |           |               */
inline TCelda<2u>::CaraVec_<TCelda<2u>::CUADR>() const          /*    |           |               */
{                                                               /*    |           |               */
TVector const u = *PtoPtrVec[1u] - *PtoPtrVec[0u],              /*    |           |               */
              v = *PtoPtrVec[2u] - *PtoPtrVec[0u];              /*    0-----------1               */

if (u[0u] * v[1u] > u[1u] * v[0u])
  return {{this, {PtoPtrVec[0u], PtoPtrVec[1u]}},
          {this, {PtoPtrVec[1u], PtoPtrVec[2u]}},
          {this, {PtoPtrVec[2u], PtoPtrVec[3u]}},
          {this, {PtoPtrVec[3u], PtoPtrVec[0u]}}};
else
  return {{this, {PtoPtrVec[0u], PtoPtrVec[3u]}},
          {this, {PtoPtrVec[3u], PtoPtrVec[2u]}},
          {this, {PtoPtrVec[2u], PtoPtrVec[1u]}},
          {this, {PtoPtrVec[1u], PtoPtrVec[0u]}}};
}

// =================================================================================================

                                                                /*               2                */
                                                                /*             ,/|`\              */
template<>                                                      /*           ,/  |  `\            */
template<>                                                      /*         ,/    '.   `\          */
std::vector<TCara<3u>>                                          /*       ,/       |     `\        */
inline TCelda<3u>::CaraVec_<TCelda<3u>::TETRA>() const          /*     ,/         |       `\      */
{                                                               /*    0-----------'.--------1     */
return {{this, {PtoPtrVec[0u], PtoPtrVec[2u], PtoPtrVec[1u]}},  /*     `\.         |      ,/      */
        {this, {PtoPtrVec[1u], PtoPtrVec[2u], PtoPtrVec[3u]}},  /*        `\.      |    ,/        */
        {this, {PtoPtrVec[0u], PtoPtrVec[3u], PtoPtrVec[2u]}},  /*           `\.   '. ,/          */
        {this, {PtoPtrVec[0u], PtoPtrVec[1u], PtoPtrVec[3u]}}}; /*              `\. |/            */
}                                                               /*                 `3             */

// =================================================================================================

template<>                                                                    /* 3----------2     */
template<>                                                                    /* |\         |\    */
std::vector<TCara<3u>>                                                        /* | \        | \   */
inline TCelda<3u>::CaraVec_<TCelda<3u>::HEXA>() const                         /* |  \       |  \  */
{                                                                             /* |   7------+---6 */
return {{this, {PtoPtrVec[3u], PtoPtrVec[2u], PtoPtrVec[1u], PtoPtrVec[0u]}}, /* |   |      |   | */
        {this, {PtoPtrVec[4u], PtoPtrVec[5u], PtoPtrVec[6u], PtoPtrVec[7u]}}, /* 0---+------1   | */
        {this, {PtoPtrVec[1u], PtoPtrVec[2u], PtoPtrVec[6u], PtoPtrVec[5u]}}, /*  \  |       \  | */
        {this, {PtoPtrVec[3u], PtoPtrVec[0u], PtoPtrVec[4u], PtoPtrVec[7u]}}, /*   \ |        \ | */
        {this, {PtoPtrVec[0u], PtoPtrVec[1u], PtoPtrVec[5u], PtoPtrVec[4u]}}, /*    \|         \| */
        {this, {PtoPtrVec[2u], PtoPtrVec[3u], PtoPtrVec[7u], PtoPtrVec[6u]}}}; /*    4----------5 */
}

// =================================================================================================

                                                                               /*        3        */
template<>                                                                     /*      ,/|`\      */
template<>                                                                     /*    ,/  |  `\    */
std::vector<TCara<3u>>                                                         /*  ,/    |    `\  */
inline TCelda<3u>::CaraVec_<TCelda<3u>::CUÑA>() const                          /* 4------+------5 */
{                                                                              /* |      |      | */
return {{this, {PtoPtrVec[2u], PtoPtrVec[1u], PtoPtrVec[0u]}},                 /* |      |      | */
        {this, {PtoPtrVec[3u], PtoPtrVec[4u], PtoPtrVec[5u]}},                 /* |      0      | */
        {this, {PtoPtrVec[0u], PtoPtrVec[1u], PtoPtrVec[4u], PtoPtrVec[3u]}},  /* |    ,/ `\    | */
        {this, {PtoPtrVec[1u], PtoPtrVec[2u], PtoPtrVec[5u], PtoPtrVec[4u]}},  /* |  ,/     `\  | */
        {this, {PtoPtrVec[0u], PtoPtrVec[3u], PtoPtrVec[5u], PtoPtrVec[2u]}}}; /* |,/         `\| */
}                                                                              /* 1-------------2 */

// =================================================================================================

                                                                /*                4               */
                                                                /*              ,/|\              */
                                                                /*            ,/ .'|\             */
                                                                /*          ,/   | | \            */
template<>                                                      /*        ,/    .' |  \           */
template<>                                                      /*      ,/      |  '.  \          */
std::vector<TCara<3u>>                                          /*    ,/        |    |  \         */
inline TCelda<3u>::CaraVec_<TCelda<3u>::PIRAM>() const          /*   0---------.'----3   \        */
{                                                               /*   `\        |      `\  \       */
return {{this, {PtoPtrVec[0u], PtoPtrVec[1u], PtoPtrVec[4u]}},  /*    `\     .'        `\  \      */
        {this, {PtoPtrVec[1u], PtoPtrVec[2u], PtoPtrVec[4u]}},  /*      `\   |           `\ \     */
        {this, {PtoPtrVec[2u], PtoPtrVec[3u], PtoPtrVec[4u]}},  /*        `\.'             `\\    */
        {this, {PtoPtrVec[4u], PtoPtrVec[3u], PtoPtrVec[0u]}},  /*           1----------------2   */
        {this, {PtoPtrVec[3u], PtoPtrVec[2u], PtoPtrVec[1u], PtoPtrVec[0u]}}};
}

// =================================================================================================

template<>
std::vector<TCara<2u>>
inline TCelda<2u>::CaraVec_() const
{
switch (Tipo)
  {
  case TRIAN : return CaraVec_<TRIAN>();
  case CUADR : return CaraVec_<CUADR>();
  default    : std::unreachable();
  }
}

// =================================================================================================

template<>
std::vector<TCara<3u>>
inline TCelda<3u>::CaraVec_() const
{
switch (Tipo)
  {
  case TETRA : return CaraVec_<TETRA>();
  case HEXA  : return CaraVec_<HEXA>();
  case CUÑA  : return CaraVec_<CUÑA>();
  case PIRAM : return CaraVec_<PIRAM>();
  default    : std::unreachable();
  }
}

// =================================================================================================

template<std::size_t d>
double
TCelda<d>::V_() const
{
double V = 0.0;

for (auto &Cara : CaraVec)
  V += Cara.Cf & Cara.Sf;
return V / d;
}

// =================================================================================================

template<std::size_t d>
TVector<d>
TCelda<d>::C_() const
{
TVector<d> C = {};

for (auto &Cara : CaraVec)
  C += Cara.Cf * (Cara.Cf & Cara.Sf);
return C / ((d + 1.0) * V);
}

// =================================================================================================
// ========================================================================================== TMalla

template<std::size_t d>
void
TMalla<d>::ReadNodes(THelper &Helper)
{
std::vector<std::size_t> NodeTagVec;
std::vector<double> CoordVec, _;
std::size_t MaxNodeTag;

gmsh::model::mesh::getNodes(NodeTagVec, CoordVec, _);
PuntoVec.reserve(NodeTagVec.size());
Helper.DefNPunto(NodeTagVec.size());
gmsh::model::mesh::getMaxNodeTag(MaxNodeTag);
Helper.DefNTag(MaxNodeTag);
for (auto const &[Tag, Coord] : std::views::zip(NodeTagVec, CoordVec | std::views::chunk(3u)))
  {
  TPunto<d> const &Punto = DefPunto(Coord);

  Helper.DefTag(Punto.ID, Tag);
  }
}

// =================================================================================================

template<std::size_t d>
template<int Dim>
void
TMalla<d>::ReadPhysGrp(std::vector<std::string> &StrVec, std::map<int, std::size_t> &TagIDMap)
{
std::vector<std::pair<int, int>> DimTagVec;

gmsh::model::getPhysicalGroups(DimTagVec, Dim);
StrVec.reserve(DimTagVec.size());
for (auto const PhysTag : DimTagVec | std::views::values)
  {
  std::string PhysName;

  gmsh::model::getPhysicalName(Dim, PhysTag, PhysName);
  TagIDMap[PhysTag] = StrVec.size();
  StrVec.emplace_back(std::move(PhysName));
  }
}

// =================================================================================================

template<std::size_t d>
template<std::size_t... i>
void
TMalla<d>::ReadElements(std::index_sequence<i...>, THelper &Helper)
{
std::map<int, std::size_t> GrpTagIDMap;
std::vector<int> ElemTypeVec;
std::vector<std::vector<std::size_t>> ElemTagVecVec,
                                      NodeTagVecVec;
std::vector<std::pair<int, int>> DimTagVec;

ReadPhysGrp<d>(GrpStrVec, GrpTagIDMap);
gmsh::model::mesh::getElements(ElemTypeVec, ElemTagVecVec, NodeTagVecVec, d);
CeldaVec.reserve(std::ranges::distance(ElemTagVecVec | std::views::join));
gmsh::model::getEntities(DimTagVec, d);
for (auto const EntityTag : DimTagVec | std::views::values)
  {
  std::vector<int> PhysTagVec;
  std::size_t GrpID = MaxID;

  gmsh::model::getPhysicalGroupsForEntity(d, EntityTag, PhysTagVec);
  if (!PhysTagVec.empty())
    GrpID = GrpTagIDMap.at(PhysTagVec.front());
  gmsh::model::mesh::getElements(ElemTypeVec, ElemTagVecVec, NodeTagVecVec, d, EntityTag);
  for (auto const &[ElemType, NodeTagVec] : std::views::zip(ElemTypeVec, NodeTagVecVec))
    {
    typename TCelda<d>::ETipo Tipo;
    std::size_t NPunto;

    switch (ElemType)
      {
      case 2  : Tipo = TCelda<d>::TRIAN; NPunto = 3u; break;
      case 3  : Tipo = TCelda<d>::CUADR; NPunto = 4u; break;
      case 4  : Tipo = TCelda<d>::TETRA; NPunto = 4u; break;
      case 5  : Tipo = TCelda<d>::HEXA;  NPunto = 8u; break;
      case 6  : Tipo = TCelda<d>::CUÑA;  NPunto = 6u; break;
      case 7  : Tipo = TCelda<d>::PIRAM; NPunto = 5u; break;
      default : continue;
      }
    for (auto const &ElemNodeTag : NodeTagVec | std::views::chunk(NPunto))
      {
      auto PtoPtrVec = ElemNodeTag
                     | std::views::transform(std::bind_front(&THelper::PuntoPtr, &Helper))
                     | std::ranges::to<std::vector>();
      TCelda<d> &Celda = DefCelda(Tipo, GrpID, std::move(PtoPtrVec));

      for (auto const Tag : ElemNodeTag)
        Helper.DefCeldaPtr(Tag, &Celda);
      }
    }
  }
for (auto &Celda : CeldaVec)
  for (auto &Cara : Celda)
    for (auto const CeldaNPtr : Helper.CeldaPtrVec(Cara.Punto(0u).ID))
      if (CeldaNPtr != &Celda && Helper.EsCeldaN(CeldaNPtr, Cara.Punto(i).ID...))
        {
        Cara.CeldaNPtr = CeldaNPtr;
        break;
        }
}

// =================================================================================================

template<std::size_t d>
template<std::size_t... i>
void
TMalla<d>::ReadBoundary(std::index_sequence<i...>, THelper const &Helper)
{
std::vector<std::pair<int, int>> DimTagVec;
std::map<int, std::size_t> CCTagIDMap;

ReadPhysGrp<d - 1>(CCStrVec, CCTagIDMap);
gmsh::model::getEntities(DimTagVec, d - 1);
for (auto const EntityTag : DimTagVec | std::views::values)
  {
  std::vector<int> PhysTagVec,
                   ElemTypeVec;
  std::vector<std::vector<std::size_t>> ElemTagVecVec,
                                        NodeTagVecVec;

  gmsh::model::getPhysicalGroupsForEntity(d - 1, EntityTag, PhysTagVec);
  if (PhysTagVec.empty())
    continue;

  std::size_t const CCID = CCTagIDMap.at(PhysTagVec.front());

  gmsh::model::mesh::getElements(ElemTypeVec, ElemTagVecVec, NodeTagVecVec, d - 1, EntityTag);
  for (auto const &[ElemType, NodeTagVec] : std::views::zip(ElemTypeVec, NodeTagVecVec))
    {
    std::size_t const NPunto = ElemType + 1u;

    for (auto const &ElemNodeTag : NodeTagVec | std::views::chunk(NPunto))
      {
      auto const PtoPtrVec = ElemNodeTag
                           | std::views::transform(std::bind_front(&THelper::PuntoPtr, &Helper))
                           | std::ranges::to<std::vector>();

      for (auto const CeldaPtr : Helper.CeldaPtrVec(PtoPtrVec.front()->ID))
        for (auto &Cara : *CeldaPtr)
          if ((... && std::ranges::contains(PtoPtrVec, &Cara.Punto(i))))
            {
            Cara.CCID = CCID;
            break;
            }
      }
    }
  }
}

// =================================================================================================

template<std::size_t d>
template<std::size_t... i>
void
TMalla<d>::ReadPeriodics(std::index_sequence<i...>, THelper const &Helper)
{
std::vector<std::pair<int, int>> DimTagVec;

gmsh::model::getEntities(DimTagVec, d - 1);
for (auto const EntityTag : DimTagVec | std::views::values)
  {
  int MasterTag;
  std::vector<std::size_t> NodeMstTagVec,
                           NodeSlvTagVec;
  std::vector<double> _;

  gmsh::model::mesh::getPeriodicNodes(d - 1, EntityTag, MasterTag, NodeMstTagVec, NodeSlvTagVec, _);
  if (EntityTag == MasterTag)
    continue;

  auto const [Beg, End] = std::ranges::subrange(NodeMstTagVec);

  for (auto It1 = Beg; It1 != End; ++It1)
    for (auto const CeldaMstPtr : Helper.CeldaPtrVec(Helper.ID(*It1)))
      for (auto &CaraMst : *CeldaMstPtr)
        {
        std::vector<std::size_t> IDVec;

        for (auto const ID : CaraMst.PtoPtrVec | std::views::transform(&TPunto<d>::ID))
          {
          auto const It2 = std::ranges::find(It1, End, ID, std::bind_front(&THelper::ID, &Helper));

          if (It2 == End)
            break;
          IDVec.push_back(Helper.ID(NodeSlvTagVec[It2 - Beg]));
          }
        if (IDVec.size() < CaraMst.NPunto())
          continue;
        for (auto const CeldaSlvPtr : Helper.CeldaPtrVec(IDVec[0u]))
          if (Helper.EsCeldaN(CeldaSlvPtr, IDVec[i]...))
            {
            for (auto &CaraSlv : *CeldaSlvPtr)
              if ((... && std::ranges::contains(CaraSlv.PtoPtrVec, IDVec[i], &TPunto<d>::ID)))
                {
                CaraSlv.CeldaNPtr = CeldaMstPtr;
                break;
                }
            CaraMst.CeldaNPtr = CeldaSlvPtr;
            break;
            }
        }
  }
}

// =================================================================================================

template<std::size_t d>
std::size_t
TMalla<d>::ID(std::vector<std::string> const &StrVec, std::string_view const Str)
{
if (auto const It = std::ranges::find(StrVec, Str); It != StrVec.end())
  return It - StrVec.begin();
throw std::invalid_argument(std::format("El grupo \"{}\" no existe", Str));
}

// =================================================================================================

template<std::size_t d>
void
TMalla<d>::Read(std::string const &FileName)
{
THelper Helper;

gmsh::initialize();
gmsh::option::setNumber("General.Terminal", 0);
gmsh::option::setNumber("Mesh.IgnorePeriodicity", 0);
gmsh::open(FileName);
ReadNodes(Helper);
ReadElements(std::make_index_sequence<d>{}, Helper);
ReadBoundary(std::make_index_sequence<d>{}, Helper);
ReadPeriodics(std::make_index_sequence<d>{}, Helper);
gmsh::finalize();

IDVecVec.resize(NCelda());
for (auto &&[IDVec, Celda] : std::views::zip(IDVecVec, CeldaVec))
  {
  std::ranges::fill(IDVec, MaxID);
  for (auto &&[ID, Cara] : std::views::zip(IDVec, Celda))
    if (!Cara.EsCC())
      {
      if (Cara.CeldaNPtr == nullptr)
        throw std::runtime_error("Malla incoherente");
      ID = Cara.CeldaN().ID;
      }
  }
}

// =================================================================================================

template<std::size_t d>
void
TMalla<d>::Write(std::ostream &os)
{
enum
{
  VTK_TRIANGLE   = 5,
  VTK_QUAD       = 9,
  VTK_TETRA      = 10,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE      = 13,
  VTK_PYRAMID    = 14
};

std::size_t CellListSize = NCelda();

os << "# vtk DataFile Version 3.0" << std::endl;
os << "VF" << std::endl;
os << "ASCII" << std::endl;
os << "DATASET UNSTRUCTURED_GRID" << std::endl;
os << "POINTS " << NPunto() << " double" << std::endl;
for (auto const &Punto : PuntoVec)
  if constexpr (d == 2u)
    os << Punto[0u] << ' ' << Punto[1u] << " 0" << std::endl;
  else
    os << Punto;
for (auto const &Celda : CeldaVec)
  CellListSize += Celda.NPunto();
os << "CELLS " << NCelda() << ' ' << CellListSize << std::endl;
for (auto const &Celda : CeldaVec)
  {
  os << Celda.NPunto();
  for (std::size_t i = 0u; i < Celda.NPunto(); ++i)
    os << ' ' << Celda.Punto(i).ID;
  os << std::endl;
  }
os << "CELL_TYPES " << NCelda() << std::endl;
for (auto const &Celda : CeldaVec)
  switch (Celda.Tipo)
    {
    case TCelda<d>::TRIAN : os << VTK_TRIANGLE   << std::endl; break;
    case TCelda<d>::CUADR : os << VTK_QUAD       << std::endl; break;
    case TCelda<d>::TETRA : os << VTK_TETRA      << std::endl; break;
    case TCelda<d>::HEXA  : os << VTK_HEXAHEDRON << std::endl; break;
    case TCelda<d>::CUÑA  : os << VTK_WEDGE      << std::endl; break;
    case TCelda<d>::PIRAM : os << VTK_PYRAMID    << std::endl;
    }
os << "CELL_DATA " << NCelda() << std::endl;
}

} // VF
