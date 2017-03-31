/****************************************************************************
  FileName  [ partitioner.cpp ]
  Synopsis  [ Implementation of the F-M two way partitioner. ]
  Author    [ Fu-Yu Chuang ]
  Date      [ 2017.3.30 ]
****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::partition()
{
    this->genInitPartition();
    this->FMAlgorithm();
    return;
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportBList()
{
    for (size_t i = 0; i < 2; ++i) {
        cout << "================ BList " << ((i == 0)? "A": "B") << "================" << endl;
        for (int j = _maxPinNum; j >= -_maxPinNum; --j) {
            cout << setw(4) << "[" << j << "] ";
            Node* node = _bList[i][j]->getNext();
            while (node != NULL) {
                cout << setw(8) << _cellArray[node->getId()]->getName() << " -> ";
                node = node->getNext();
            }
            cout << endl;
        }
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

// Private member functions
void Partitioner::genInitPartition()
{
    bool part = 0;
    int count = 1, mode = 7;
    int tmpNet = _cellArray[0]->getFirstNet();
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        if (_cellArray[i]->getFirstNet() == tmpNet) {
            _cellArray[i]->setPart(part);
            ++_partSize[part];
        }
        else {
            if (count % mode != 0) {
                tmpNet = _cellArray[i]->getFirstNet();
                _cellArray[i]->setPart(part);
                ++_partSize[part];
                ++count;
            }
            else {
                part = !part;
                tmpNet = _cellArray[i]->getFirstNet();
                _cellArray[i]->setPart(part);
                ++_partSize[part];
                count = 1;
            }
        }
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j)
            _netArray[netList[j]]->incPartCount(part);
    }

    // Make sure the initial partiiton is balanced
    while (!this->checkBalance()) {
        this->reBalance();
    }

    this->initGain();
    this->buildBList();
    this->countCutSize();
    return;
}

void Partitioner::FMAlgorithm()
{
    bool run = true;
    this->initPass();

    while (run && (_moveNum < _cellNum)) {
        if (_unlockNum[0] == 0) {
            Cell* max = this->findMaxGainCell(1);
            if (abs(_partSize[1] - _partSize[0] - 2) < (_bFactor * _cellNum))
                updateGain(max);
            else
                run = false;

        }
        else if (_unlockNum[1] == 0) {
            Cell* max = this->findMaxGainCell(0);
            if (abs(_partSize[0] - _partSize[1] - 2) < (_bFactor * _cellNum))
                updateGain(max);
            else
                run = false;
        }
        else {
            Cell* maxA = this->findMaxGainCell(0);
            Cell* maxB = this->findMaxGainCell(1);
            if (maxA->getGain() >= maxB->getGain()) {
                if (abs(_partSize[0] - _partSize[1] - 2) < (_bFactor * _cellNum))
                    updateGain(maxA);
                else if (abs(_partSize[1] - _partSize[0] - 2) < (_bFactor * _cellNum))
                    updateGain(maxB);
                else
                    run = false;
            }
            else {
                if (abs(_partSize[1] - _partSize[0] - 2) < (_bFactor * _cellNum))
                    updateGain(maxB);
                else if (abs(_partSize[0] - _partSize[1] - 2) < (_bFactor * _cellNum))
                    updateGain(maxA);
                else
                    run = false;
            }
        }
        ++_moveNum;
    }

    if (_maxAccGain > 0) {
        ++_iterNum;

        this->recover2Best();
        cout << endl;
        cout << "------------------------------------------" << endl;
        cout << " Pass #" << _iterNum << endl;
        cout << " Gain: " << _maxAccGain << endl;
        cout << "------------------------------------------" << endl;

        this->FMAlgorithm();
    }
    this->countCutSize();
    return;
}

void Partitioner::insertCell(Cell* c)
{
    int gain = c->getGain();
    bool part = c->getPart();
    Node* node = c->getNode();
    node->setPrev(_bList[part][gain]);
    node->setNext(_bList[part][gain]->getNext());
    _bList[part][gain]->setNext(node);
    if (node->getNext() != NULL)
        node->getNext()->setPrev(node);
    return;
}

void Partitioner::removeCell(Cell* c)
{
    Node* node = c->getNode();
    node->getPrev()->setNext(node->getNext());
    if (node->getNext() != NULL)
        node->getNext()->setPrev(node->getPrev());
    return;
}

// called when gain is updated
void Partitioner::moveCell(Cell* c)
{
    removeCell(c);
    insertCell(c);
    return;
}

void Partitioner::buildBList()
{
    _bList[0].clear();
    _bList[1].clear();

    this->countMaxPinNum();
    for (int i = -_maxPinNum; i <= _maxPinNum; ++i) {
        if (_bList[0].count(i) == 0)
            _bList[0][i] = new Node(-1);    // dummy node
        if (_bList[1].count(i) == 0)
            _bList[1][i] = new Node(-1);    // dummy node
    }

    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        this->insertCell(_cellArray[i]);
    }
    return;
}

void Partitioner::initGain()
{
    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];

    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        bool part = _cellArray[i]->getPart();
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            if (_netArray[netList[j]]->getPartCount(int(part)) == 1)
                _cellArray[i]->incGain();
            if (_netArray[netList[j]]->getPartCount(int(!part)) == 0)
                _cellArray[i]->decGain();
        }
    }
    return;
}

void Partitioner::updateGain(Cell* c)
{
    _accGain += c->getGain();

    bool fPart = c->getPart();
    bool tPart = !c->getPart();
    c->lock();
    _moveStack.push_back(c->getNode()->getId());

    vector<int> netList = c->getNetList();
    for (size_t i = 0, end_i = netList.size(); i < end_i; ++i) {
        if (_netArray[netList[i]]->getPartCount(tPart) == 0) {
            vector<int> cellList = _netArray[netList[i]]->getCellList();
            for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
                if (!_cellArray[cellList[j]]->getLock()) {
                    _cellArray[cellList[j]]->incGain();
                    this->moveCell(_cellArray[cellList[j]]);
                }
            }
        }
        else if (_netArray[netList[i]]->getPartCount(tPart) == 1) {
            vector<int> cellList = _netArray[netList[i]]->getCellList();
            for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
                if (!_cellArray[cellList[j]]->getLock() &&
                    _cellArray[cellList[j]]->getPart() == tPart) {
                    _cellArray[cellList[j]]->decGain();
                    this->moveCell(_cellArray[cellList[j]]);
                }
            }
        }

        _netArray[netList[i]]->decPartCount(fPart);
        _netArray[netList[i]]->incPartCount(tPart);
        c->setPart(tPart);

        if (_netArray[netList[i]]->getPartCount(fPart) == 0) {
            vector<int> cellList = _netArray[netList[i]]->getCellList();
            for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
                if (!_cellArray[cellList[j]]->getLock()) {
                    _cellArray[cellList[j]]->decGain();
                    this->moveCell(_cellArray[cellList[j]]);
                }
            }
        }
        else if (_netArray[netList[i]]->getPartCount(fPart) == 1) {
            vector<int> cellList = _netArray[netList[i]]->getCellList();
            for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
                if (!_cellArray[cellList[j]]->getLock() &&
                    _cellArray[cellList[j]]->getPart() == fPart) {
                    _cellArray[cellList[j]]->incGain();
                    this->moveCell(_cellArray[cellList[j]]);
                }
            }
        }
    }
    this->removeCell(c);
    ++_partSize[tPart];
    --_partSize[fPart];
    --_unlockNum[fPart];
    if (_accGain > _maxAccGain) {
        this->storeBestState();
    }
    return;
}

Cell* Partitioner::findMaxGainCell(bool part)
{
    int maxGain = _maxPinNum;
    while (maxGain >= -_maxPinNum && _bList[part][maxGain]->getNext() == NULL) {
        --maxGain;
    }
    Cell* maxGainCell = _cellArray[_bList[part][maxGain]->getNext()->getId()];

    return maxGainCell;
}

void Partitioner::initPass()
{
    for (size_t i = 0; i < _cellNum; ++i) {
        _cellArray[i]->unlock();
        _cellArray[i]->setGain(0);
    }

    this->initGain();
    this->buildBList();
    _accGain = 0;
    _maxAccGain = 0;
    _moveNum = 0;
    _bestMoveNum = 0;
    _moveStack.clear();
    return;
}

void Partitioner::countCutSize()
{
    _cutSize = 0;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        if (_netArray[i]->getPartCount(0) && _netArray[i]->getPartCount(1)) {
            ++_cutSize;
        }
    }
    return;
}

void Partitioner::countMaxPinNum()
{
    _maxPinNum = 0;
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPinNum() > _maxPinNum) {
            _maxPinNum = _cellArray[i]->getPinNum();
        }
    }
    return;
}

void Partitioner::reBalance()
{
    cout << "Re-balancing the initial partition..." << endl;
    int diff = abs(_partSize[0]-_partSize[1]);
    bool biggerPart = (_partSize[0] > _partSize[1])? 0: 1;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        if (_cellArray[i]->getPart() == biggerPart) {
            _cellArray[i]->setPart(!biggerPart);
            diff -= 2;
            ++_partSize[!biggerPart];
            --_partSize[biggerPart];
            vector<int> netList = _cellArray[i]->getNetList();
            for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
                _netArray[netList[j]]->incPartCount(!biggerPart);
                _netArray[netList[j]]->decPartCount(biggerPart);
            }
            if (diff <= 0)
                break;
        }
    }
    return;
}

bool Partitioner::checkBalance()
{
    bool balanced = true;
    if (_partSize[0] < ((1 - _bFactor) / 2) * _cellNum ||
        _partSize[0] > ((1 + _bFactor) / 2) * _cellNum ||
        _partSize[1] < ((1 - _bFactor) / 2) * _cellNum ||
        _partSize[1] > ((1 + _bFactor) / 2) * _cellNum) {
        balanced = false;
    }
    return balanced;
}

void Partitioner::storeBestState()
{
    _maxAccGain = _accGain;
    _bestMoveNum = _moveNum;
    return;
}

void Partitioner::recover2Best()
{
    _moveNum = _bestMoveNum;
    for (size_t i = _moveStack.size()-1; i > _moveNum; --i) {
        bool part = _cellArray[_moveStack[i]]->getPart();
        _cellArray[_moveStack[i]]->setPart(!part);
        --_partSize[part];
        ++_partSize[!part];
    }
    for (size_t i = 0; i < _netNum; ++i) {
        _netArray[i]->setPartCount(0, 0);
        _netArray[i]->setPartCount(1, 0);
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end = cellList.size(); j < end; ++j) {
            _netArray[i]->incPartCount(_cellArray[cellList[j]]->getPart());
        }
    }
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
