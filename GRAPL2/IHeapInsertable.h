#pragma once
template<typename KeyType> class IHeapInsertable
{
public:
	IHeapInsertable<KeyType>()
	{
		key = new KeyType;
	}

	virtual ~IHeapInsertable()
	{
		delete key;
	}

	virtual bool operator <(const IHeapInsertable<KeyType>& other) const
	{
		return *key < *(other.key);
	}

	virtual bool operator <=(const IHeapInsertable<KeyType>& other) const
	{
		return *key <= *(other.key);
	}

	virtual bool operator >(const IHeapInsertable<KeyType>& other) const
	{
		return *key > *(other.key);
	}

	virtual bool operator >=(const IHeapInsertable<KeyType>& other) const
	{
		return *key >= *(other.key);
	}

	virtual bool operator ==(const IHeapInsertable<KeyType>& other) const
	{
		return *key == *(other.key);
	}


	virtual void OnIndexChanged(unsigned int NewIndex) = 0;
	virtual void ModifyKey(KeyType NewKey) = 0;

	KeyType* key;
};