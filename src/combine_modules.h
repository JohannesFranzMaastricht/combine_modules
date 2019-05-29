// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef COMBINE_MODULES_H_
#define COMBINE_MODULES_H_

#include "biodynamo.h"
#include "biology-modules.h"
#include "core/substance_initializers.h"
#include "my-neurite.h"
#include "my-cell.h"
#include "neuroscience/neuroscience.h"

namespace bdm {

// Define compile time parameter
BDM_CTPARAM(experimental::neuroscience) {
  BDM_CTPARAM_HEADER(experimental::neuroscience);

  using NeuronSoma = MyCell;
  using NeuriteElement = MyNeurite;

  using SimObjectTypes = CTList<bdm::MyCell, bdm::MyNeurite>;

  BDM_CTPARAM_FOR(bdm, MyCell) { 
    //using BiologyModules = CTList<NullBiologyModule>;
    using BiologyModules = CTList<GrowDivide>; 
  };

  BDM_CTPARAM_FOR(bdm, MyNeurite) {
    using BiologyModules = CTList<ApicalElongationBM, BasalElongationBM>;
  };
};

inline int Simulate(int argc, const char** argv) {
  Simulation<> simulation(argc, argv);

  auto* rm = simulation.GetResourceManager();
  rm->template Reserve<MyCell>(1);

  auto construct_soma = [](const std::array<double, 3> position) {
    MyCell cell(position);
    cell.SetDiameter(6);
    cell.SetAdherence(0.4);
    cell.SetMass(1.0);
    cell.SetCanDivide(true);
    cell.AddBiologyModule(GrowDivide(8, 400, {gAllEventIds})); // 
    return cell;
  };

  auto add_dendrites = [](MyCell soma) {
    auto soma_soptr = soma.GetSoPtr();
    auto&& dendrite_apical = soma_soptr->ExtendNewNeurite({0, 0, 1});
    dendrite_apical->AddBiologyModule(ApicalElongationBM());
    dendrite_apical->SetCanBranch(true);
    auto&& dendrite_basal1 = soma_soptr->ExtendNewNeurite({0, 0, -1});
    dendrite_basal1->AddBiologyModule(BasalElongationBM());
  };

  auto cell = construct_soma({0, 0, 0});
  rm->push_back(cell);
  add_dendrites(cell);

  // Run simulation for one timestep
  simulation.GetScheduler()->Simulate(10);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace bdm

#endif  // COMBINE_MODULES_H_
